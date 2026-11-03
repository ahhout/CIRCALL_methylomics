#!/usr/bin/env python3
"""Verify shipped UXM pooling and prepare real CpG counts for a matched CELFIE panel."""
import argparse
import concurrent.futures
import csv
import datetime as dt
import gzip
import itertools
import json
import os
from pathlib import Path
import sys

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from output_contract import output_guard

from pat_counts import atomic_json, compile_counter, count_pat, sha, utc

ATLAS_SHA = 'e58fbec4fe6200a311540c6e1a16be901a2d1e9ed9d74b6a57888f2feb635f51'
DONOR_AUDIT = 'docs/audits/2026-10-26_12-33-uxm-reference-donor-audit-v2/specimen_inventory.tsv'
MK_CANDIDATE = 'docs/audits/2026-10-26_12-58-megakaryocyte-count-investigation-v1/reference_candidate.json'


def read_tsv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def load_counts(path):
    records = read_tsv(path)
    return np.array([[int(row[key]) for key in ('methylated', 'total', 'U', 'X', 'M')] for row in records], dtype=np.int64)


def check_pool(atlas_rows, cell, counts):
    expected = np.array([float(row[cell]) if row[cell] != 'NA' else np.nan for row in atlas_rows])
    depth = counts[:, 2:5].sum(axis=1)
    numerator = np.where([row['direction'] == 'U' for row in atlas_rows], counts[:, 2], counts[:, 4])
    score = np.divide(numerator, depth, out=np.full(len(depth), np.nan), where=depth > 0)
    reconstructed = np.round(score, 3)
    matches = (reconstructed == expected) | (np.isnan(reconstructed) & np.isnan(expected))
    finite = np.isfinite(expected)
    return {'all_values_match': bool(matches.all()), 'finite_atlas_values': int(finite.sum()),
            'finite_matches': int((matches & finite).sum()), 'missing_atlas_values': int((~finite).sum()),
            'missing_matches': int((matches & ~finite).sum()),
            'mismatch_rows_1based': (np.where(~matches)[0] + 1).tolist()}


def write_matrix(path, atlas_rows, names, arrays):
    stage = path.with_name(path.name + '.partial')
    with stage.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['chr', 'start', 'end'] + [name + suffix for name in names for suffix in ('_METH', '_DEPTH')])
        for index, region in enumerate(atlas_rows):
            writer.writerow([region[key] for key in ('chr', 'start', 'end')] +
                            [int(array[index, col]) for array in arrays for col in (0, 1)])
    stage.replace(path)


def run(args):
    project = args.project.resolve()
    revision = output_guard.require_clean_source(Path(__file__))
    output_guard.require_clean_source(project)
    atlas = project / 'workflows/uxm/reference/Atlas.U25.l4.hg38.full.tsv'
    if sha(atlas) != ATLAS_SHA:
        raise ValueError('The existing l4/40-cell atlas checksum differs')
    atlas_rows = read_tsv(atlas)
    cells = list(atlas_rows[0])[8:]
    if len(atlas_rows) != 979 or len(cells) != 40:
        raise ValueError('Expected unchanged atlas with 979 regions and 40 cells')
    output = args.output.absolute()
    dt.datetime.strptime(output.name[:16], '%Y-%m-%d_%H-%M')
    if len(output.name) < 18 or output.name[16] != '-':
        raise ValueError('Output name must start YYYY-MM-DD_HH-MM-')
    metadata = project / (DONOR_AUDIT if args.phase == 'reference' else 'workflows/deconflow/config/samples.pats.tsv')
    records = read_tsv(metadata)
    if args.phase == 'reference':
        inputs = [{'sample': row['gsm_id'], 'pat': str(project / 'references/uxm_atlas_pat' / row['pat_filename']),
                   'expected_size': int(row['pat_size_bytes'])} for row in records]
    else:
        inputs = [{'sample': row['sample'], 'pat': row['pat']} for row in records]
    if len({item['sample'] for item in inputs}) != len(inputs):
        raise ValueError('Duplicate count input IDs')
    for item in inputs:
        pat = Path(item['pat'])
        if not pat.is_file():
            raise ValueError(f'Missing PAT {pat}')
        if 'expected_size' in item and pat.stat().st_size != item['expected_size']:
            raise ValueError(f'PAT size does not match recorded inventory: {pat}')
    settings = {'project_revision': revision, 'phase': args.phase, 'atlas_sha256': ATLAS_SHA, 'atlas': str(atlas),
                'input_manifest': str(metadata), 'input_manifest_sha256': sha(metadata),
                'helper_sha256': sha(__file__), 'counter_helper_sha256': sha(Path(__file__).with_name('pat_counts.py')),
                'coordinate_convention': 'atlas coordinates unchanged; CpG membership [startCpG,endCpG)',
                'count_semantics': 'actual integer C/T CpG observations, multiplicity weighted; no minimum fragment length',
                'reference_pool_verification': 'native wgbstools homog strict l4 U/X/M; pooled then rounded to three decimals',
                'inputs': inputs}
    execution = output / 'execution.json'
    if getattr(args, 'resume', False):
        output_guard.resume(output, settings)
        previous = json.loads(execution.read_text())
        if previous['settings'] != settings or previous['status'] in ('completed', 'blocked_pool_mismatch'):
            raise ValueError('Cannot resume changed or completed count preparation')
    else:
        output_guard.reserve_path(output, settings)
        previous = None
    output_guard.update(output, status='running')
    provenance = previous or {'started_utc': utc(), 'settings': settings, 'attempts': []}
    provenance['status'] = 'running'
    provenance['attempts'].append({'started_utc': utc(), 'workers': args.workers})
    atomic_json(execution, provenance)
    compiled = compile_counter(output / 'counter')
    try:
        def counter(item):
            result = count_pat(item['pat'], atlas, output / 'per_sample' / item['sample'], args.tools_dir, compiled=compiled)
            print(json.dumps({'sample': item['sample'], 'status': result['status'], 'stream': result['stream']}), flush=True)
            return item['sample'], result
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            counts_provenance = dict(executor.map(counter, inputs))
        arrays = {item['sample']: load_counts(output / 'per_sample' / item['sample'] / 'counts.tsv') for item in inputs}
        if args.phase == 'patients':
            target = output / 'patient.counts.tsv'
            names = [item['sample'] for item in inputs]
            write_matrix(target, atlas_rows, names, [arrays[name] for name in names])
            gate = {'status': 'passed', 'atlas_sha256': ATLAS_SHA, 'patient_counts_sha256': sha(target),
                    'patient_counts': str(target), 'sample_order': names, 'regions': len(atlas_rows),
                    'started_utc': provenance['started_utc'], 'completed_utc': utc()}
        else:
            grouping = {}
            for record in records:
                cell = record['candidate_atlas_cell']
                if cell:
                    grouping.setdefault(cell, []).append(record['gsm_id'])
            # Full-stream all-six heart checks replace the earlier padded indexed investigation.
            grouping['Heart-Cardio'] = sorted([item['sample'] for item in inputs if item['sample'] in {
                'GSM5652212','GSM5652213','GSM5652214','GSM5652215','GSM5652216','GSM5652217'}])
            validation = {}
            pools = {}
            selected = {}
            for cell in cells:
                if cell == 'Megakaryocytes':
                    continue
                members = grouping.get(cell, [])
                if not members:
                    raise ValueError(f'No PAT candidate group for {cell}')
                full = sum((arrays[name] for name in members), np.zeros((len(atlas_rows), 5), dtype=np.int64))
                result = check_pool(atlas_rows, cell, full)
                matching = []
                examined = 1
                if cell == 'Heart-Cardio' or (cell == 'Ovary-Ep' and not result['all_values_match']):
                    examined = 0
                    for length in range(1, len(members) + 1):
                        for subset in itertools.combinations(members, length):
                            candidate = sum((arrays[name] for name in subset), np.zeros_like(full))
                            check = check_pool(atlas_rows, cell, candidate)
                            examined += 1
                            if check['all_values_match']:
                                matching.append((list(subset), candidate, check))
                    if len(matching) == 1:
                        members, full, result = matching[0]
                    elif len(matching) > 1:
                        result['all_values_match'] = False
                        result['ambiguity'] = 'Multiple numerically identical candidate subsets; no pool selected'
                result.update(selected_samples=members, candidate_subsets_examined=examined,
                              matching_subsets=[value[0] for value in matching],
                              evidence='Exact equality of native count-derived rounded scores including missingness; historical uniqueness only within tested candidate pools')
                validation[cell] = result
                pools[cell] = full
                selected[cell] = members
                print(json.dumps({'cell': cell, 'validation': result}), flush=True)
            atomic_json(output / 'pool_validation.json', validation)
            mk_manifest_path = project / MK_CANDIDATE
            mk = json.loads(mk_manifest_path.read_text())
            mk_path = Path(mk['counts']['path'])
            if sha(mk_path) != mk['counts']['sha256']:
                raise ValueError('Prior MK component checksum differs')
            donor_counts = {str(donor['donor_number']): np.zeros((len(atlas_rows), 2), dtype=np.int64) for donor in mk['donor_profiles']}
            with gzip.open(mk_path, 'rt') as handle:
                for record in csv.DictReader(handle, delimiter='\t'):
                    if record['candidate_pool'] == mk['selected_candidate']:
                        index = int(record['atlas_row_1based']) - 1
                        donor_counts[record['donor_number']][index] = [int(record['methylated']), int(record['total'])]
            pools['Megakaryocytes'] = sum(donor_counts.values(), np.zeros((len(atlas_rows), 2), dtype=np.int64))
            atomic_json(output / 'reference_pooling.json', {
                'pat_cell_pools': selected,
                'megakaryocyte_component': str(mk_manifest_path), 'megakaryocyte_manifest_sha256': sha(mk_manifest_path),
                'megakaryocyte_counts_sha256': sha(mk_path),
                'megakaryocyte_policy': 'Sum actual counts across documented all71-library donor profiles; preserves donor profiles in original component; original UXM fragment pooling remains unestablished',
                'missing_coverage': 'Keep every region and cell; METH=DEPTH=0 explicitly represents no observed CpGs'})
            all_pass = all(value['all_values_match'] for value in validation.values())
            gate = {'status': 'passed' if all_pass else 'blocked_pool_mismatch', 'atlas_sha256': ATLAS_SHA,
                    'pool_validation_sha256': sha(output / 'pool_validation.json'), 'pat_cell_columns_checked': len(validation),
                    'pat_cell_columns_passed': sum(value['all_values_match'] for value in validation.values()),
                    'unresolved_cells': [cell for cell, value in validation.items() if not value['all_values_match']],
                    'megakaryocyte_original_uxm_pattern_equivalence': 'not established; documented count reconstruction',
                    'started_utc': provenance['started_utc'], 'completed_utc': utc()}
            if all_pass:
                target = output / 'reference.counts.tsv'
                write_matrix(target, atlas_rows, cells, [pools[cell] for cell in cells])
                gate.update(reference_counts_sha256=sha(target), reference_counts=str(target), cell_order=cells,
                            regions=len(atlas_rows), zero_depth_regions_by_cell={cell: int((pools[cell][:, 1] == 0).sum()) for cell in cells})
        atomic_json(output / 'validation_gate.json', gate)
        provenance.update(status='completed' if gate['status'] == 'passed' else 'blocked_pool_mismatch', completed_utc=utc(),
                          validation_gate=gate, counts_provenance={name: {'compressed_pat_sha256': value['compressed_pat_sha256'],
                                                                      'counts_sha256': value['counts_sha256']} for name, value in counts_provenance.items()})
        provenance['attempts'][-1]['completed_utc'] = utc()
        atomic_json(execution, provenance)
        output_guard.finish(output, True)
        print(json.dumps(gate, indent=2), flush=True)
        return 0 if gate['status'] == 'passed' else 2
    except BaseException as error:
        provenance['status'] = 'failed'
        provenance['attempts'][-1].update(failed_utc=utc(), error=repr(error))
        atomic_json(execution, provenance)
        output_guard.finish(output, False)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--phase', choices=('reference', 'patients'), required=True)
    parser.add_argument('--resume', action='store_true', help='Continue only this owned incomplete preparation')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--tools-dir', type=Path, default=Path('/path/to/hpc-home/.conda/envs/nextflow'))
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Full count preparation must run through sbatch in a Slurm allocation')
    if args.workers < 1:
        parser.error('--workers must be positive')
    sys.exit(run(args))


if __name__ == '__main__':
    main()
