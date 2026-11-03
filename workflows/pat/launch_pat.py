#!/usr/bin/env python3
"""Submit-side PAT bundle orchestration; run through run_prepare_pat.sh in Slurm."""
import argparse
import csv
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import signal
import shutil
import subprocess
import sys

from prepare_pat import sha256
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from output_contract import output_guard, identity

REPO = Path(__file__).resolve().parents[2]
SOURCE = REPO / 'workflows/pat'
RESULTS = Path('/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/uxm_deconvolution')
REFERENCE = Path('/path/to/hpc-home/.conda/envs/nextflow/share/wgbs_tools')


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('results-dir', 'metadata', 'samplesheet'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--run-name', required=True, help='New base name; output appends -pat-v1')
    p.add_argument('--conversion-mode', choices=('reuse', 'complete'), required=True)
    p.add_argument('--reuse-dir', type=Path, action='append', default=[])
    p.add_argument('--reuse-sample', action='append', default=[],
                   help='In complete mode, explicitly select samples to reuse; others are generated')
    return p.parse_args(argv)


def samples(path, delimiter):
    with path.open(newline='') as h:
        rows = list(csv.DictReader(h, delimiter=delimiter))
    ids = [r.get('sample', '') for r in rows]
    if not ids or any(not x for x in ids) or len(ids) != len(set(ids)):
        raise ValueError(f'Nonempty unique sample IDs required: {path}')
    return ids


def absent(path):
    if path.exists() or path.is_symlink():
        raise ValueError(f'Output/work destination exists; inspect without overwriting: {path}')


def command(argv, log=None):
    if log:
        with log.open('x') as h:
            subprocess.run([str(x) for x in argv], check=True, stdout=h)
    else:
        subprocess.run([str(x) for x in argv], check=True)


def launch(args):
    started_utc = os.environ.get('PAT_EXECUTION_STARTED_UTC') or datetime.now(timezone.utc).isoformat()
    args._execution_path = None
    args._owned_output = None
    try:
        _launch(args, started_utc)
    except BaseException:
        if args._owned_output is not None:
            output_guard.finish(args._owned_output, False)
        if args._execution_path is not None:
            try:
                record = json.loads(args._execution_path.read_text())
                record.pop('completed_utc', None)
                record.update(status='failed_or_interrupted', ended_utc=datetime.now(timezone.utc).isoformat())
                args._execution_path.write_text(json.dumps(record, indent=2) + '\n')
            except (OSError, ValueError):
                pass  # Retain the original error if failure provenance cannot be written.
        raise


def _launch(args, started_utc):
    job = os.environ.get('SLURM_JOB_ID', '')
    if not job.isdigit():
        raise ValueError('Submit with sbatch; numeric SLURM_JOB_ID required')
    match = re.fullmatch(r'(\d{4}-\d{2}-\d{2}_[0-2][0-9]-[0-5][0-9])-[A-Za-z0-9][A-Za-z0-9._-]*', args.run_name)
    if not match:
        raise ValueError('Run name must start with YYYY-MM-DD_HH-MM (UTC), followed by a description')
    date_text = match.group(1)
    datetime.strptime(date_text, '%Y-%m-%d_%H-%M')
    threads = int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))
    if threads < 1:
        raise ValueError('Positive CPU allocation required')
    args.results_dir = args.results_dir.resolve(strict=True)
    args.metadata = args.metadata.resolve(strict=True)
    args.samplesheet = args.samplesheet.resolve(strict=True)
    args.reuse_dir = [p.resolve(strict=True) for p in args.reuse_dir]
    ids = samples(args.metadata, '\t')
    if ids != samples(args.samplesheet, ','):
        raise ValueError('Metadata and samplesheet must contain identical samples in identical order')
    if len(args.reuse_sample) != len(set(args.reuse_sample)) or not set(args.reuse_sample) <= set(ids):
        raise ValueError('Reuse samples must be unique selected sample IDs')
    if args.conversion_mode == 'reuse' and (not args.reuse_dir or args.reuse_sample):
        raise ValueError('Reuse mode requires bundles and reuses the whole cohort; omit --reuse-sample')
    if args.conversion_mode == 'complete' and bool(args.reuse_dir) != bool(args.reuse_sample):
        raise ValueError('Complete mode requires both reuse bundles and reuse samples, or neither')
    final = RESULTS / (args.run_name + '-pat-v1')
    work = RESULTS / f'.{args.run_name}-pat-v1.running-{job}'
    absent(final); absent(work)
    required = [SOURCE / n for n in ('prepare_pat.py', 'validate_pat_manifest.py', 'launch_pat.py', 'run_prepare_pat.sh')]
    required += [REPO / n for n in ('environment.yml', 'config/nextflow.config', 'config/params.emseq.integrated.yml', 'samplesheets/samplesheet.csv', 'references/emseq_hg38_lambda_puc19/manifest.tsv')]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'Missing required file: {path}')
    if not (args.results_dir / 'bismark/deduplicated').is_dir():
        raise ValueError('Missing source bismark/deduplicated directory')
    executables = {}
    for name in ('wgbstools', 'samtools', 'tabix'):
        exe = shutil.which(name)
        if not exe:
            raise ValueError(f'Missing executable: {name}')
        executables[name] = exe
    if not RESULTS.is_dir() or not os.access(RESULTS, os.W_OK):
        raise ValueError(f'Missing or unwritable publication root: {RESULTS}')
    # Stable lock persists; it is never unlinked while another process may hold it.
    with (RESULTS / (final.name + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        absent(final); absent(work)
        output_guard.reserve_path(final, identity(SOURCE, [args.metadata, args.samplesheet], work))
        args._owned_output = final
        work.mkdir()
        args._execution_path = work / 'publication_execution.json'
        args._execution_path.write_text(json.dumps({
            'started_utc': started_utc, 'status': 'running', 'job_id': job}, indent=2) + '\n')
        submitted = work / 'submitted'; submitted.mkdir()
        snapshots = required + [args.metadata, args.samplesheet]
        destinations = [p.name for p in required] + ['metadata.tsv', 'input_samplesheet.csv']
        for src, name in zip(snapshots, destinations):
            shutil.copyfile(src, submitted / name)
        (submitted / 'invocation.json').write_text(json.dumps({
            'job_id': job, 'source_results': str(args.results_dir),
            'metadata': str(args.metadata), 'samplesheet': str(args.samplesheet),
            'run_name': args.run_name, 'conversion_mode': args.conversion_mode,
            'reuse_dirs': [str(p) for p in args.reuse_dir], 'reuse_samples': args.reuse_sample,
            'threads': threads, 'source_revision': subprocess.check_output(
                ['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip()
        }, indent=2) + '\n')
        (work / 'submitted.sha256').write_text(''.join(
            f'{sha256(p)}  submitted/{p.name}\n' for p in sorted(submitted.iterdir())))
        bundle = work / 'bundle'
        prepare = [sys.executable, submitted / 'prepare_pat.py', '--results-dir', args.results_dir,
                   '--metadata', submitted / 'metadata.tsv', '--output-dir', bundle, '--owner-output', final,
                   '--conversion-mode', args.conversion_mode, '--genome', 'hg38',
                   '--threads', str(threads), '--reference-root', REFERENCE,
                   '--reference-manifest', REPO / 'references/emseq_hg38_lambda_puc19/manifest.tsv']
        for name, exe in executables.items(): prepare += ['--' + name, exe]
        for directory in args.reuse_dir: prepare += ['--reuse-dir', directory]
        for sample in args.reuse_sample: prepare += ['--reuse-sample', sample]
        command(prepare)
        validate = [sys.executable, submitted / 'validate_pat_manifest.py', '--results-dir', args.results_dir,
                    '--samplesheet', submitted / 'input_samplesheet.csv', '--reference-root', REFERENCE,
                    '--tabix', executables['tabix']]
        command(validate + ['--pat-dir', bundle], work / 'validation_before_publish.log')
        for name in ('submitted', 'submitted.sha256', 'validation_before_publish.log', 'publication_execution.json'):
            (work / name).rename(bundle / name)
        args._execution_path = bundle / 'publication_execution.json'
        # Move individual artifacts without replacing anything in the owned destination.
        output_guard.publish_contents(bundle, final)
        args._execution_path = final / 'publication_execution.json'
        validate[1] = final / 'submitted/validate_pat_manifest.py'
        validate[validate.index('--samplesheet') + 1] = final / 'submitted/input_samplesheet.csv'
        command(validate + ['--pat-dir', final], final / 'validation_after_publish.log')
        record = json.loads(args._execution_path.read_text())
        ended_utc = datetime.now(timezone.utc).isoformat()
        record.update(status='completed', completed_utc=ended_utc, ended_utc=ended_utc)
        args._execution_path.write_text(json.dumps(record, indent=2) + '\n')
        output_guard.finish(final, True)
        print(f'PAT bundle published and validated: {final}')


def main():
    args = parse_args()
    try:
        launch(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f'ERROR: {exc}') from exc


if __name__ == '__main__':
    def interrupted(signum, frame):
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, interrupted)
    main()
