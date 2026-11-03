#!/usr/bin/env python3
"""Validate frozen inputs and record provenance; no methylation percentages used."""
import argparse
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import shutil


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    tmp = Path(str(path) + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2) + '\n')
    os.replace(tmp, path)


def validate_run_name(path):
    if not re.fullmatch(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-[A-Za-z0-9][A-Za-z0-9._-]*', path.name):
        raise ValueError('Expected date-first UTC run name YYYY-MM-DD_HH-MM-purpose-vN')
    dt.datetime.strptime(path.name[:16], '%Y-%m-%d_%H-%M')


def inference_sha(path):
    return sha(path)


def validate_samples(rows):
    if len(rows) != 12 or len({r['sample'] for r in rows}) != 12:
        raise ValueError('Expected twelve unique original sample IDs')
    if len({str(Path(r['report']).resolve()) for r in rows}) != 12:
        raise ValueError('Each sample requires a distinct source report')
    groups = collections.defaultdict(list)
    for r in rows:
        if not re.fullmatch(r'[A-Za-z0-9.-]+', r['sample']):
            raise ValueError('Unsafe sample ID')
        if Path(r['report']).name != r['sample'] + '.CpG_report.txt.gz':
            raise ValueError('Sample/report mismatch')
        if r['upstream_run'] not in Path(r['report']).parts:
            raise ValueError('Upstream run/report mismatch')
        if r['reference_name'] != 'wgbstools_hg38_lambda_puc19_v1':
            raise ValueError('Reference mismatch')
        groups[r['patient']].append(r['compartment'])
    if len(groups) != 4 or any(sorted(v) != ['BM', 'PB', 'cfDNA'] for v in groups.values()):
        raise ValueError('Expected four complete, unique patient triplets')


def scan_report(path, interval, destination):
    totals = collections.defaultdict(lambda: [0, 0, 0])
    selected = collections.defaultdict(lambda: [0, 0])
    out = gzip.open(destination, 'wt') if interval else None
    try:
        with gzip.open(path, 'rt') as src:
            for line in src:
                x = line.rstrip('\n').split('\t')
                if len(x) != 7 or x[2] not in ('+', '-') or x[5] != 'CG':
                    raise ValueError(f'Invalid Bismark CpG record in {path}')
                pos, m, u = int(x[1]), int(x[3]), int(x[4])
                if pos < 1 or min(m, u) < 0 or max(m, u, m + u) > 2147483647:
                    raise ValueError('Invalid integer count/coordinate')
                canonical = pos - (x[2] == '-')
                if canonical < 1:
                    raise ValueError('Reverse-strand coordinate outside reference')
                t = totals[x[0]]; t[0] += 1; t[1] += m; t[2] += m + u
                if interval and x[0] == interval[0] and interval[1] <= canonical <= interval[2]:
                    out.write(line)
                    k = (x[0], canonical); selected[k][0] += m; selected[k][1] += m + u
    finally:
        if out:
            out.close()
    return totals, selected


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', required=True, type=Path)
    p.add_argument('--run', required=True, type=Path)
    p.add_argument('--resume', action='store_true')
    a = p.parse_args(); c = json.loads(a.config.read_text()); run = a.run
    validate_run_name(run)
    if run.is_symlink() or not run.is_dir():
        raise ValueError('Create the run through run_dmrcate.sh first')
    if (run / 'COMPLETE').exists():
        raise ValueError('Refusing completed bundle')
    import sys
    project = Path('/path/to/hpc-home/analysis/CIRCALL_methylomics')
    sys.path.insert(0, str(project / 'scripts/lib'))
    import output_guard
    revision = output_guard.require_clean_source(project)
    if sha(__file__) != sha(project / 'workflows/dmrcate/prepare.py'):
        raise ValueError('Submitted preparation helper differs from committed source')
    ownership = output_guard.read_record(run)
    if ownership['identity'].get('project_revision') != revision:
        raise ValueError('Preparation requires the original project revision')
    if ownership.get('output') != str(run.resolve()) or ownership.get('status') not in ('running', 'failed'):
        raise ValueError('Preparation requires an owned incomplete run from the new launcher')
    with (run / 'preparation_attempts.jsonl').open('a') as log:
        log.write(json.dumps({'started_utc': utc(), 'resume': a.resume}) + '\n')
    with open(c['samples']) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    validate_samples(rows)
    for path in c['upstream_params']:
        upstream = json.loads(Path(path).read_text())
        if (Path(upstream['fasta']).resolve() != Path(c['reference_fasta']).resolve()
                or upstream['aligner'] != 'bismark' or upstream['em_seq'] is not True):
            raise ValueError('Upstream alignment/reference provenance mismatch')
    runtime = Path(c['runtime'])
    if not (runtime / 'STAGED').is_file():
        raise ValueError('Runtime and annotation staging is incomplete')
    expected = (runtime / 'sif.sha256').read_text().split()[0]
    if sha(runtime / 'dmrcate.sif') != expected:
        raise ValueError('SIF checksum mismatch')
    if sha(c['reference_fasta']) != c['reference_sha256']:
        raise ValueError('Reference sequence differs from approved provenance')
    with open(c['reference_manifest']) as f:
        reference = {r['field']: r['value'] for r in csv.DictReader(f, delimiter='\t')}
    if (reference['combined_fasta_sha256'] != c['reference_sha256'] or
            sha(reference['source_controls']) != reference['source_controls_sha256']):
        raise ValueError('Reference/control manifest checksum mismatch')
    subprocess.run(['sha256sum', '-c', 'SHA256SUMS'], cwd=runtime / 'annotations', check=True)
    sources = [Path(__file__), Path(__file__).with_name('analysis.R'), a.config, Path(c['samples']),
               Path(c['reference_manifest']), Path(c['namespace']), Path(c['chrom_sizes']),
               runtime / 'packages.lock.tsv', runtime / 'annotations/SHA256SUMS']
    sources += [Path(x) for x in c['upstream_params'] + c['control_qc']]
    sources += [Path(reference['source_controls'])]
    fingerprints = {str(x.resolve()): sha(x) for x in sources}
    fingerprints['runtime_sif'] = expected
    fingerprints['inference_code'] = inference_sha(Path(__file__).with_name('analysis.R'))
    for r in rows:
        path = Path(r['report'])
        if not path.is_file():
            raise ValueError(f'Missing report: {path}')
        fingerprints[str(path)] = sha(path)
    marker = run / 'input_fingerprint.json'
    if marker.exists():
        if not a.resume or json.loads(marker.read_text()) != fingerprints:
            raise ValueError('Resume requires identical inputs, configuration, runtime and inference code')
    else:
        write_json(marker, fingerprints)
    if (run / 'prepared/READY').exists():
        with (run / 'preparation_attempts.jsonl').open('a') as log:
            log.write(json.dumps({'completed_utc': utc(), 'reused': True}) + '\n')
        return
    prep = run / 'prepared'; prep.mkdir(exist_ok=True)
    provenance = run / 'provenance'; provenance.mkdir(exist_ok=True)
    for i, path in enumerate(c['control_qc']):
        shutil.copyfile(path, provenance / f'upstream_control_qc_{i+1}.tsv')
    for i, path in enumerate(c['upstream_params']):
        shutil.copyfile(path, provenance / f'upstream_params_{i+1}.json')
    for path in (c['samples'], c['reference_manifest'], runtime / 'packages.lock.tsv',
                 runtime / 'sources/downloads.tsv', runtime / 'annotations/sources.tsv'):
        label = 'annotation_sources.tsv' if str(path).endswith('annotations/sources.tsv') else Path(path).name
        shutil.copyfile(path, provenance / label)
    # Incomplete attempts remain as recovery evidence; a new preparation directory
    # is required after failure during import (no silent overwrite of HDF5 files).
    control_rows = []; check_rows = []
    for r in rows:
        output = prep / (r['sample'] + '.CpG_report.txt.gz')
        totals, selected = scan_report(r['report'], c['smoke_interval'], output)
        r['analysis_report'] = str(output) if c['smoke_interval'] else r['report']
        for chrom, values in totals.items():
            control_rows.append([r['sample'], chrom] + values)
        for (chrom, pos), (m, cov) in selected.items():
            check_rows.append([r['sample'], chrom, pos, m, cov])
    with (prep / 'samples.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
    for name, header, data in [('source_count_totals.tsv', ['sample','chr','strand_records','M','Cov'], control_rows),
                               ('smoke_expected_counts.tsv', ['sample','chr','pos','M','Cov'], check_rows)]:
        with (prep / name).open('w') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(header); w.writerows(data)
    write_json(prep / 'READY', {'completed_utc': utc()})
    ownership = json.loads((run / '.run-ownership.json').read_text())
    if ownership.get('output') != str(run.resolve()) or ownership.get('status') not in ('running', 'failed'):
        raise ValueError('Preparation requires an owned incomplete run from the new launcher')
    with (run / 'preparation_attempts.jsonl').open('a') as log:
        log.write(json.dumps({'completed_utc': utc(), 'reused': False}) + '\n')


if __name__ == '__main__':
    main()
