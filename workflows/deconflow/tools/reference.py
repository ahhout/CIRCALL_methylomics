#!/usr/bin/env python3
"""Audit published membership and build a separate beta reference only after it passes."""
import argparse
import os
import sys
import xml.etree.ElementTree as ET
import zipfile
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import *


def validate_marker_coordinates(data):
    """Check both interval endpoints against the pinned one-based CpG dictionary."""
    positions = {(r['chr'], int(r[key]) + offset)
                 for r in data for key, offset in [('start', 0), ('end', -1)]}
    result = run(['tabix', GENOME / 'CpG.bed.gz'] +
                 [f'{chrom}:{pos}-{pos}' for chrom, pos in sorted(positions)],
                 capture_output=True, text=True)
    lookup = {(chrom, int(pos)): int(index)
              for chrom, pos, index in (line.split() for line in result.stdout.splitlines())}
    for marker in data:
        if (int(marker['startCpG']) >= int(marker['endCpG']) or
            lookup.get((marker['chr'], int(marker['start']))) != int(marker['startCpG']) or
            lookup.get((marker['chr'], int(marker['end']) - 1)) != int(marker['endCpG']) - 1):
            raise ValueError(f'Marker coordinates disagree with the canonical CpG dictionary: {marker["name"]}')
    return len(data)


def published_samples(workbook):
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(workbook) as z:
        strings = [''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si', ns)]
        records = []
        for row in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:row', ns)[3:]:
            data = {}
            for c in row.findall('s:c', ns):
                value = c.find('s:v', ns)
                value = value.text if value is not None else ''
                if c.get('t') == 's':
                    value = strings[int(value)]
                data[re.sub('[0-9]', '', c.get('r'))] = value
            if data.get('A'):
                records.append({'published_sample': data['A'], 'sample_key': data['A'].rsplit('-', 1)[-1],
                                'published_group': data['E'], 'refined_group': data['F'],
                                'source_row': row.get('r')})
    if len(records) != 205 or len({r['sample_key'] for r in records}) != 205:
        raise ValueError('Unexpected Table S1 shape or nonunique specimen identifiers')
    return records


def audit(root):
    source = read_json(WORKFLOW / 'config/source.lock.json')
    workbook = WORKFLOW / 'software/loyfer_supplement.xlsx'
    if sha(workbook) != source['supplement_sha256']:
        raise ValueError('Published workbook changed')
    _, cells = atlas_rows()
    published = {r['sample_key']: r for r in published_samples(workbook)}
    manifest = PROJECT / 'references/uxm_atlas_pat/download_manifest.tsv'
    local = rows(manifest)
    result = []
    for entry in local:
        p = published.get(entry['donor_id'])
        group = p['published_group'] if p else ''
        # Exact published group matches only. Abbreviations/composite labels are
        # explicitly left for review instead of assuming atlas training membership.
        status = 'verified' if group in cells else 'unresolved'
        pat = manifest.parent / entry['pat_filename']
        csi = manifest.parent / entry['csi_filename']
        if not pat.is_file() or not csi.is_file():
            raise ValueError(f'Missing reference payload: {entry["gsm_id"]}')
        if pat.stat().st_size != int(entry['pat_size_bytes']) or csi.stat().st_size != int(entry['csi_size_bytes']):
            raise ValueError(f'Reference size mismatch: {entry["gsm_id"]}')
        result.append(dict(entry, **(p or {}), atlas_cell_type=group if status == 'verified' else '',
                           status=status, pat=str(pat), csi=str(csi), evidence_url=source['supplement_url']))
    columns = ['gsm_id', 'sample_name', 'donor_id', 'published_sample', 'published_group', 'refined_group',
               'atlas_cell_type', 'status', 'pat', 'csi', 'source_row', 'evidence_url']
    table(root / 'membership.tsv', result, columns)
    observed = {r['atlas_cell_type'] for r in result if r['status'] == 'verified'}
    missing = sorted(set(cells) - observed)
    pending = [r['gsm_id'] for r in result if r['status'] != 'verified']
    report = {'started_utc': read_json(root / 'execution.json')['started_utc'], 'checked_utc': utc(),
              'atlas': identity(ATLAS, True), 'supplement': identity(workbook, True),
              'download_manifest': identity(manifest, True), 'membership_sha256': sha(root / 'membership.tsv'),
              'required_cells': cells, 'missing_cells': missing, 'unresolved_samples': pending,
              'published_samples': 205, 'downloaded_samples': len(local),
              'status': 'blocked' if missing or pending else 'passed',
              'note': 'GSE186458 lacks megakaryocytes. No cell types dropped, aliases assumed, samples excluded or U scores converted.'}
    write_json(root / 'reference_audit.json', report)
    if missing or pending:
        finish(root, 'blocked', missing_cells=missing, unresolved_samples=pending)
        return 2
    finish(root)
    (root / 'MEMBERSHIP_VALIDATED').write_text(sha(root / 'membership.tsv') + '\n')
    return 0


def membership(path, cells):
    data = rows(path)
    if not data or len({r['gsm_id'] for r in data}) != len(data):
        raise ValueError('Empty or duplicate reference accession mapping')
    if any(r['status'] not in ('verified', 'excluded-reviewed') or not r['evidence_url'] for r in data):
        raise ValueError('Unresolved reference membership; review evidence before building')
    selected = [r for r in data if r['status'] == 'verified']
    if set(r['atlas_cell_type'] for r in selected) != set(cells):
        raise ValueError('Reference must cover exactly the approved 40 atlas cell types')
    if len({r['pat'] for r in selected}) != len(selected):
        raise ValueError('Reference PAT assigned more than once')
    return selected


def equal_sample_beta(counts, groups, cells):
    """counts is samples × markers × (methylated,total); missing is not zero."""
    counts = np.asarray(counts)
    if counts.ndim != 3 or counts.shape[2] != 2 or len(groups) != counts.shape[0]:
        raise ValueError('Invalid regional count dimensions')
    if not np.isfinite(counts).all() or (counts < 0).any() or (counts[..., 0] > counts[..., 1]).any():
        raise ValueError('Invalid methylated/total counts')
    beta = np.divide(counts[..., 0], counts[..., 1], out=np.full(counts.shape[:2], np.nan), where=counts[..., 1] > 0)
    matrix, donors = [], []
    for cell in cells:
        values = beta[np.array(groups) == cell]
        n = np.isfinite(values).sum(axis=0)
        if (n == 0).any():
            raise ValueError(f'Missing reference beta for {cell}; no imputation permitted')
        matrix.append(np.nansum(values, axis=0) / n)
        donors.append(n)
    return np.array(matrix).T, np.array(donors).T


def build(root, mapping):
    data, cells = atlas_rows()
    selected = membership(mapping, cells)  # gate BEFORE any count conversion
    version = run(['wgbstools', '--version'], capture_output=True, text=True).stdout.strip()
    if version != 'wgbstools version 0.3.0':
        raise ValueError('Reference preparation requires existing wgbstools 0.3.0')
    if sha(GENOME / 'CpG.bed.gz') != CPG_SHA:
        raise ValueError('Reference namespace mismatch')
    validate_marker_coordinates(data)
    tools_path = PROJECT.parent / '.conda/envs/nextflow/share/wgbs_tools/src/python'
    sys.path.insert(0, str(tools_path))
    from beta_to_blocks import collapse_process, is_block_file_nice
    blocks = pd.DataFrame(data)[['chr', 'start', 'end', 'startCpG', 'endCpG']]
    for col in blocks.columns[1:]:
        blocks[col] = blocks[col].astype(int)
    blocks.to_csv(root / 'markers.tsv', sep='\t', index=False)
    nice, _ = is_block_file_nice(blocks)
    (root / 'counts').mkdir()
    all_counts, provenance = [], []
    for record in selected:
        pat, csi = Path(record['pat']), Path(record['csi'])
        run(['tabix', '-l', pat], root / 'commands.jsonl', capture_output=True)
        out = root / 'counts' / record['gsm_id']
        out.mkdir()
        run(['wgbstools', 'pat2beta', pat, '--genome', 'hg38', '--lbeta', '-@', '4', '-o', out], root / 'commands.jsonl')
        beta = out / (pat.name.removesuffix('.pat.gz') + '.lbeta')
        array = np.memmap(beta, dtype=np.uint16, mode='r')
        if len(array) != 29401795 * 2 or (array == 65535).any():
            raise ValueError('Wrong namespace size or saturated LBETA counts; stop without changing scientific behavior')
        count = collapse_process(str(beta), blocks, nice)
        if count is None:
            raise ValueError('wgbstools regional aggregation failed')
        all_counts.append(count)
        pd.DataFrame(count, columns=['methylated', 'total']).to_csv(out / 'regional_counts.tsv', sep='\t', index=False)
        provenance.append(dict(record, pat_identity=identity(pat, True), csi_identity=identity(csi, True), lbeta=identity(beta, True)))
    values, n = equal_sample_beta(np.array(all_counts), [r['atlas_cell_type'] for r in selected], cells)
    for filename, array in [('reference.beta.tsv', values), ('reference.contributors.tsv', n)]:
        frame = pd.concat([blocks.iloc[:, :3], pd.DataFrame(array, columns=cells)], axis=1)
        frame.to_csv(root / filename, sep='\t', index=False, float_format='%.12g')
    write_json(root / 'reference.provenance.json', {'atlas': identity(ATLAS, True), 'membership': identity(mapping, True),
               'inputs': provenance, 'aggregation': 'regional M/total per sample; equal sample means per cell type',
               'minimum_observations': 1, 'imputation': False, 'dictionary_sha256': CPG_SHA,
               'wgbstools_version': version,
               'wgbstools_source': {name: sha(tools_path / name) for name in ['pat2beta.py', 'beta_to_blocks.py', 'utils_wgbs.py']},
               'beta_sha256': sha(root / 'reference.beta.tsv')})
    finish(root)
    (root / 'REFERENCE_VALIDATED').write_text(sha(root / 'reference.beta.tsv') + '\n')
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['audit', 'build'])
    parser.add_argument('output', type=Path)
    parser.add_argument('--membership', type=Path)
    args = parser.parse_args()
    if args.action == 'build' and not os.environ.get('SLURM_JOB_ID'):
        parser.error('Reference conversion must run in a Slurm allocation')
    if args.action == 'build' and not args.membership:
        parser.error('--membership is required')
    root = new_run(args.output)
    try:
        return audit(root) if args.action == 'audit' else build(root, args.membership)
    except Exception as exc:
        finish(root, 'failed', error=str(exc))
        raise


if __name__ == '__main__':
    sys.exit(main())
