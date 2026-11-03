#!/usr/bin/env python3
"""Freeze EpiDISH recovery inputs and validate stage completion receipts."""
import argparse
import csv
import hashlib
import json
import os
import platform
from importlib.metadata import version
import subprocess
from datetime import datetime, timezone
from pathlib import Path

PREPARED_FILES = (
    'marker_mask.tsv', 'coverage_summary.tsv', 'sample_regional_values.tsv.gz',
    'reference_matrix.tsv', 'epidish_rpc_fractions.tsv', 'epidish_warnings.tsv',
    'uxm_native_fractions.tsv', 'uxm_matched_atlas.l4.tsv',
)

FINAL_COMPARISON_FILES = (
    'uxm_weighted_matched_fractions.tsv', 'comparison_metrics.tsv',
    'composition_comparison.pdf', 'fraction_scatter.pdf', 'manifest.json',
)
FINAL_SUMMARY_FILES = (
    'fractions_long.tsv.gz', 'comparison_metrics_all.tsv', 'epidish_warnings_all.tsv',
    'threshold_stability.tsv', 'threshold_stability.pdf', 'marker_counts.tsv', 'analysis_report.md',
)
FINAL_EPIDISH_FILES = (
    'epidish_full.pdf', 'epidish_full.png', 'epidish_hematopoietic.pdf', 'epidish_hematopoietic.png',
    'epidish_annotated_fractions.tsv', 'epidish_plot_provenance.json', 'fractions_long.tsv.gz',
    'epidish_warnings_all.tsv', 'marker_counts.tsv',
)


def finalization_outputs():
    return ([f'comparisons/{profile}/coverage-ge{threshold}/{name}'
             for profile in ('full', 'hematopoietic') for threshold in (1, 3, 5)
             for name in FINAL_COMPARISON_FILES] +
            [f'comparisons/summary/{name}' for name in FINAL_SUMMARY_FILES] +
            [f'epidish_r/{name}' for name in FINAL_EPIDISH_FILES] +
            ['provenance/sessionInfo.txt', 'provenance/run_manifest.json'])


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds').replace('+00:00', 'Z')


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def record(path):
    path = Path(path).resolve(strict=True)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing or empty contract input: {path}')
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha256_file(path)}


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    os.replace(temporary, path)


def receipt_path(root, stage):
    return Path(root) / 'provenance' / 'receipts' / (stage.replace('/', '__') + '.json')


def input_fingerprint(root, stage):
    value = sha256_file(Path(root) / 'provenance' / 'recovery_contract.json')
    if stage.startswith(('matched/', 'nnls/')):
        value += ':' + sha256_file(receipt_path(root, 'prepared/' + stage.split('/', 1)[1]))
    return value


def safe_relative(root, relative):
    path = Path(relative)
    if path.is_absolute() or '..' in path.parts or not path.parts:
        raise ValueError(f'Unsafe receipt output: {relative}')
    return Path(root) / path


def check_receipt(root, stage, outputs):
    receipt = receipt_path(root, stage)
    paths = [safe_relative(root, name) for name in outputs]
    if not receipt.exists():
        if any(os.path.lexists(path) for path in paths):
            raise ValueError(f'Unreceipted partial stage {stage}; use a fresh run directory')
        return False
    data = json.loads(receipt.read_text())
    if data.get('stage') != stage or data.get('input_fingerprint') != input_fingerprint(root, stage):
        raise ValueError(f'Receipt input mismatch for {stage}')
    records = data.get('outputs', [])
    if [r.get('path') for r in records] != list(outputs):
        raise ValueError(f'Receipt output list mismatch for {stage}')
    for path, evidence in zip(paths, records):
        if not path.is_file() or path.stat().st_size != evidence['bytes'] or sha256_file(path) != evidence['sha256']:
            raise ValueError(f'Receipt output checksum mismatch: {path}')
    return True


def write_receipt(root, stage, outputs):
    target = receipt_path(root, stage)
    if target.exists():
        raise ValueError(f'Receipt already exists: {target}')
    evidence = []
    for name in outputs:
        path = safe_relative(root, name)
        evidence.append({'path': name, 'bytes': path.stat().st_size, 'sha256': sha256_file(path)})
    target.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(target, {'stage': stage, 'input_fingerprint': input_fingerprint(root, stage),
                         'completed_utc': utc_now(), 'outputs': evidence})


def check_recovery_outputs(root):
    root = Path(root)
    # All reusable outputs must have a complete, valid receipt before any resume write.
    for profile in ('full', 'hematopoietic'):
        for threshold in (1, 3, 5):
            suffix = f'{profile}/coverage-ge{threshold}'
            directory = f'comparisons/{suffix}'
            check_receipt(root, f'prepared/{suffix}', [f'{directory}/{n}' for n in PREPARED_FILES])
            check_receipt(root, f'matched/{suffix}', [f'{directory}/.uxm_weighted_matched.raw.csv'])
            check_receipt(root, f'nnls/{suffix}', [f'{directory}/nnls_unweighted_fractions.tsv'])
    preparation = ['provenance/comparison_jobs.tsv', 'provenance/pat_manifest.tsv',
                   'provenance/sample_metadata.tsv', 'provenance/input_manifest.json']
    check_receipt(root, 'preparation', preparation)
    check_receipt(root, 'nnls_provenance', ['provenance/unweighted_nnls.json'])
    check_receipt(root, 'matched_reuse_provenance', ['provenance/matched_uxm_reuse.json'])
    if check_receipt(root, 'finalization', finalization_outputs()):
        manifest = json.loads((root / 'provenance/run_manifest.json').read_text())
        linked = {'input_manifest_sha256': 'input_manifest.json',
                  'commands_log_sha256': 'commands.log',
                  'unweighted_nnls_sha256': 'unweighted_nnls.json',
                  'recovery_contract_sha256': 'recovery_contract.json'}
        if manifest.get('matched_uxm_reuse_sha256') is not None:
            linked['matched_uxm_reuse_sha256'] = 'matched_uxm_reuse.json'
        for field, name in linked.items():
            if manifest.get(field) != sha256_file(root / 'provenance' / name):
                raise ValueError(f'Finalized provenance checksum mismatch: {name}')


def build_contract(args):
    full, heme = Path(args.full_run), Path(args.heme_run)
    manifests = []
    for run in (full, heme):
        with (run / 'pat_manifest.tsv').open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        if not rows or list(rows[0]) != ['sample', 'pat', 'csi', 'lbeta']:
            raise ValueError(f'Invalid PAT manifest in {run}')
        manifests.append(rows)
    if manifests[0] != manifests[1]:
        raise ValueError('Full/heme PAT manifests differ')
    samples = [row['sample'] for row in manifests[0]]
    if len(set(samples)) != len(samples):
        raise ValueError('Duplicate PAT samples')
    artifacts = []
    for row in manifests[0]:
        for kind in ('pat', 'csi', 'lbeta'):
            artifacts.append({'sample': row['sample'], 'kind': kind, **record(safe_relative(full, row[kind]))})
    regional = []
    for sample in samples:
        candidates = list((full / 'intermediates' / 'l4').glob(f'{sample}.mem.*.homog.gz'))
        if len(candidates) != 1:
            raise ValueError(f'Expected one regional table for {sample}')
        regional.append({'sample': sample, **record(candidates[0])})
    scripts = Path(__file__).resolve().parent
    identities = [scripts.parent / 'run_epidish_comparison.sh', scripts.parent / 'environment.yml',
                  scripts.parents[2] / 'environment.yml']
    identities.extend(sorted(scripts.glob('*.R')))
    identities.extend(sorted(scripts.glob('*.py')))
    baseline = None
    if args.baseline:
        base = Path(args.baseline).resolve(strict=True)
        if not (base / 'EPIDISHR_COMPLETED').is_file():
            raise ValueError('Baseline is not completed')
        baseline = {'directory': str(base), 'run_manifest': record(base / 'provenance/run_manifest.json')}
    return {'schema': 1, 'atlas': record(args.atlas), 'metadata': record(args.metadata),
            'samplesheet': record(args.samplesheet), 'samples': samples,
            'profiles': ['full', 'hematopoietic'], 'coverage_thresholds': [1, 3, 5],
            'regional_measure': 'U/(U+X+M) for U markers; M/(U+X+M) for M markers',
            'epidish_method': 'RPC', 'epidish_maxit': args.rpc_maxit,
            'r_identity': args.r_identity, 'python_identity': platform.python_version(),
            'python_packages': {name: version(name) for name in ('numpy', 'pandas', 'scipy')},
            'uxm_commit': subprocess.check_output(['git', '-C', args.uxm_checkout, 'rev-parse', 'HEAD'], text=True).strip(),
            'full_manifest': record(full / 'pat_manifest.tsv'),
            'heme_manifest': record(heme / 'pat_manifest.tsv'),
            'native_results': [record(run / 'results/lineage_fractions.csv') for run in (full, heme)],
            'regional_tables': regional, 'pat_artifacts': artifacts,
            'script_environment_identities': [record(path) for path in identities],
            'baseline': baseline, 'source_results': str(Path(args.source_results).resolve(strict=True))}


def reserve(args):
    root = Path(args.work_dir)
    path = root / 'provenance' / 'recovery_contract.json'
    if args.recover and not path.is_file():
        raise ValueError('Missing recovery contract; legacy incomplete work requires a fresh run directory')
    contract = build_contract(args)
    if args.recover:
        if json.loads(path.read_text()) != contract:
            raise ValueError('Recovery contract differs; use the original launcher revision and inputs or a fresh run directory')
        check_recovery_outputs(root)
        timing = json.loads((root / 'provenance/execution_attempts.json').read_text())
        if not timing.get('started_utc') or not timing.get('attempts'):
            raise ValueError('Missing original execution timing; use a fresh run directory')
        if any(item.get('attempt_id') == args.attempt_id for item in timing['attempts']):
            raise ValueError('Execution attempt ID already exists; each resume needs a distinct job ID')
    else:
        root.mkdir()  # exclusive reservation: never adopt pre-existing work
        (root / 'provenance').mkdir()
        atomic_json(path, contract)
        timing = {'started_utc': args.started_utc, 'attempts': []}
    timing['attempts'].append({'attempt_id': args.attempt_id, 'started_utc': args.started_utc, 'status': 'running'})
    atomic_json(root / 'provenance/execution_attempts.json', timing)


def finish(args):
    path = Path(args.work_dir) / 'provenance/execution_attempts.json'
    data = json.loads(path.read_text())
    attempt = data['attempts'][-1]
    if attempt['attempt_id'] != args.attempt_id:
        raise ValueError('Refusing to update another execution attempt')
    attempt['status'] = args.status
    attempt['ended_utc'] = utc_now()
    if args.status == 'completed':
        attempt['completed_utc'] = attempt['ended_utc']
    else:
        attempt.pop('completed_utc', None)
    atomic_json(path, data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    freeze = sub.add_parser('reserve')
    for name in ('work-dir', 'atlas', 'metadata', 'samplesheet', 'full-run', 'heme-run', 'uxm-checkout', 'source-results', 'r-identity', 'started-utc', 'attempt-id'):
        freeze.add_argument('--' + name, required=True)
    freeze.add_argument('--rpc-maxit', type=int, required=True)
    freeze.add_argument('--baseline', default='')
    freeze.add_argument('--recover', action='store_true')
    done = sub.add_parser('finish')
    done.add_argument('--work-dir', required=True)
    done.add_argument('--attempt-id', required=True)
    done.add_argument('--status', choices=('completed', 'failed'), required=True)
    for operation in ('check', 'receipt'):
        command = sub.add_parser(operation)
        command.add_argument('work_dir')
        command.add_argument('stage')
        command.add_argument('outputs', nargs='+')
    args = parser.parse_args()
    try:
        if args.command == 'reserve':
            reserve(args)
        elif args.command == 'finish':
            finish(args)
        elif args.command == 'receipt':
            write_receipt(args.work_dir, args.stage, args.outputs)
        elif not check_receipt(args.work_dir, args.stage, args.outputs):
            return 3
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f'EpiDISH recovery refused: {error}\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
