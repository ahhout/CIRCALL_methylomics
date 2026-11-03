#!/usr/bin/env python3
"""Fit patient CelFiE from validated actual counts on the fixed hg38 UXM regions."""
import argparse
from contextlib import contextmanager
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
from common import (ATLAS, ATLAS_SHA, COMMIT, WORK, WORKFLOW, atlas_rows, finish,
                    identity, new_run, require_incomplete, output_guard, resume_nextflow_session, save_nextflow_session, read_json, rows, run, sha, utc, write_json)
from native import checked, verify_runtime

SETTINGS = dict(celfie_maxiter=1000, unknown=0, converg=0.0001,
                celfie_randrest=10, parall_job=1)
OWN_FILES = ['tools/celfie_patients.py', 'tools/celfie_patients.nf',
             'tools/celfie_diagnostics.py', 'tools/run_celfie_patients.sh',
             'common.py', 'native.py']


@contextmanager
def run_lock(root):
    with (root / '.lock').open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Another CelFiE attempt is already running in this output directory') from None
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def runtime_versions():
    commands = {'nextflow': ['nextflow', '-version'],
                'apptainer': ['apptainer', '--version']}
    return {name: subprocess.check_output(command, text=True).strip()
            for name, command in commands.items()}


def count_table(path, expected_coordinates, labels=None):
    """Validate integer counts without substituting or dropping any coordinates."""
    with Path(path).open(newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        records = list(reader)
    if len(header) != len(set(header)) or [c.lower() for c in header[:3]] != ['chr', 'start', 'end']:
        raise ValueError(f'Invalid coordinate or duplicate count header: {path}')
    if (len(header) - 3) % 2 or len(header) <= 3:
        raise ValueError(f'Expected methylated/depth column pairs: {path}')
    found = []
    for i in range(3, len(header), 2):
        methylated, depth = header[i:i + 2]
        if not methylated.lower().endswith('_meth'):
            raise ValueError(f'Invalid methylated column: {methylated}')
        name = methylated[:-5]
        if depth.lower() != (name + '_depth').lower() or not name:
            raise ValueError(f'Invalid methylated/depth pair: {methylated}, {depth}')
        found.append(name)
    if len(found) != len(set(found)) or (labels is not None and found != labels):
        raise ValueError(f'Unexpected count labels or order: {found}')
    if len(records) != len(expected_coordinates):
        raise ValueError(f'Unexpected region count: {path}')
    for row, coordinates in zip(records, expected_coordinates):
        if len(row) != len(header) or row[:3] != list(coordinates):
            raise ValueError(f'Coordinates or row width differ from fixed atlas: {path}')
        for i in range(3, len(row), 2):
            pair = row[i:i + 2]
            if any(not re.fullmatch(r'0|[1-9]\d*', value) for value in pair):
                raise ValueError(f'Counts must be nonnegative integers: {path}, {coordinates}')
            methylated, depth = map(int, pair)
            if methylated > depth or depth > 2 ** 53:
                raise ValueError(f'Invalid or inexact native count range: {path}, {coordinates}')
    frame = pd.DataFrame(records, columns=header)
    return frame, found


def native_table(reference, patients, cells, samples):
    """Native CelFiE expects a second coordinate triplet before reference pairs."""
    if not samples or len(samples) != len(set(samples)):
        raise ValueError('Supply a nonempty unique sample list')
    if set(samples) & set(cells):
        raise ValueError('Sample and reference labels overlap')
    if any('_' in cell for cell in cells):
        raise ValueError('Pinned native tissue-name parser cannot retain underscores')
    output = patients.iloc[:, :3].copy()
    output.columns = ['CHR', 'START', 'END']
    for frame, selected, prefix in ((patients, samples, ''), (reference, cells, 'REF_')):
        if prefix:
            for i, coordinate in enumerate(('CHR', 'START', 'END')):
                output[prefix + coordinate] = reference.iloc[:, i]
        for label in selected:
            for suffix in ('meth', 'depth'):
                matches = [c for c in frame.columns[3:] if c.lower() == (label + '_' + suffix).lower()]
                if len(matches) != 1:
                    raise ValueError(f'Missing or ambiguous count column: {label}_{suffix}')
                output[label + '_' + suffix] = frame[matches[0]].values
            if not any(int(value) > 0 for value in output[label + '_depth']):
                raise ValueError(f'No observed counts for selected sample/reference: {label}')
    return output.to_csv(sep='\t', index=False)


def gate_reference(path, reference, patients):
    gate = read_json(path)
    if gate.get('status') != 'passed':
        raise ValueError('Reference validation gate has not passed')
    if gate.get('atlas_sha256') != ATLAS_SHA:
        raise ValueError('Reference gate does not identify the fixed UXM atlas')
    if gate.get('reference_counts_sha256') != sha(reference):
        raise ValueError('Reference counts differ from the validated reference')
    if 'patient_counts_sha256' in gate and gate['patient_counts_sha256'] != sha(patients):
        raise ValueError('Patient counts differ from the validated patient counts')
    return gate


def immutable_text(path, text):
    path = Path(path)
    if path.exists():
        if path.read_text() != text:
            raise ValueError(f'Refusing changed run-local input/configuration: {path}')
    else:
        with path.open('x') as handle:
            handle.write(text)


def render_config(runtime, root, entrypoint):
    for path in (runtime, root, WORK, entrypoint):
        if any(c in str(path) for c in "'\n,\\"):
            raise ValueError(f'Unsupported character in Nextflow path: {path}')
    source = runtime / ('DecoNFlow-' + COMMIT)
    return f"""// Frozen run-local configuration for the unchanged native CelFiE process.
process {{
    executor = 'slurm'
    queue = 'shortq'
    cpus = 2
    memory = '4 GB'
    time = '2h'
    errorStrategy = 'terminate'
    maxRetries = 0
    container = '{runtime}/images/python.sif'
    beforeScript = 'mkdir -p .tmp'
    containerOptions = {{ '--bind ' + task.workDir + '/.tmp:/tmp --bind {runtime}:{runtime}:ro --env TMPDIR=/tmp' }}
    publishDir = [path: '{root}/results', mode: 'copy', overwrite: false]
    withName: RUN_DECONV {{ container = '{runtime}/images/python.sif' }}
    withName: OBSERVE_RESTARTS {{ publishDir = [path: '{root}/diagnostics', mode: 'copy', overwrite: false] }}
}}
executor {{ queueSize = 2; submitRateLimit = '1/1s' }}
apptainer {{ enabled = true; autoMounts = true; cacheDir = '/path/to/hpc-scratch/nextflow/apptainer_cache' }}
env {{
    PATH = '{source}/bin:/usr/local/bin:/usr/bin:/bin'
    TZ = 'UTC'
    TMPDIR = '/tmp'
    OMP_NUM_THREADS = '1'
    OPENBLAS_NUM_THREADS = '1'
    http_proxy = 'http://127.0.0.1:9'
    https_proxy = 'http://127.0.0.1:9'
    HTTP_PROXY = 'http://127.0.0.1:9'
    HTTPS_PROXY = 'http://127.0.0.1:9'
}}
trace {{
    enabled = true
    file = params.trace_path
    fields = 'task_id,hash,native_id,name,status,exit,duration,peak_rss,workdir'
}}
cleanup = false
"""


def validate_results(root, samples, cells, diagnostics):
    output = root / 'results/1_tissue_proportions.txt'
    frame = checked(pd.read_csv(output, sep='\t', index_col=0), samples, cells)
    result = dict(fractions=identity(output, digest=True),
                  samples=samples, cells=cells, regions=979,
                  minimum_fraction=float(frame.values.min()), maximum_fraction=float(frame.values.max()))
    if diagnostics:
        replay = root / 'diagnostics/1_tissue_proportions.txt'
        if sha(output) != sha(replay):
            raise ValueError('Observer replay changed the native selected fractions')
        document = read_json(root / 'diagnostics/restart_diagnostics.json')
        if document['samples'] != samples or document['cells'] != cells:
            raise ValueError('Restart diagnostics have unexpected labels')
        if len(document['restarts']) != SETTINGS['celfie_randrest']:
            raise ValueError('Incomplete native restart diagnostics')
        if not document['best_start_converged']:
            raise ValueError('Selected native CelFiE restart did not converge')
        result['restart_diagnostics'] = identity(root / 'diagnostics/restart_diagnostics.json', digest=True)
        result['all_starts_converged'] = document['all_starts_converged']
        result['best_start_converged'] = document['best_start_converged']
        result['observer_fractions_identical'] = True
    return result


def invoke(root, runtime, params, phase, attempt, resume=False, require_cache=False):
    require_incomplete(root)
    trace = root / f'attempt-{attempt:03d}-{phase}.trace.tsv'
    command = ['nextflow', '-C', root / 'celfie.config', 'run', WORKFLOW / 'tools/celfie_patients.nf',
               '-params-file', root / 'params.json',
               '-work-dir', WORK, '--trace_path', trace]
    if resume:
        command.extend(['-resume', resume_nextflow_session(root)])
    env = dict(os.environ, NXF_OFFLINE='true', NXF_ANSI_LOG='false',
               NXF_ENABLE_VIRTUAL_THREADS='false',
               NXF_PLUGINS_DIR=str(runtime / 'plugins'), NXF_WORK=str(WORK))
    with (root / f'attempt-{attempt:03d}-{phase}.log').open('x') as handle:
        try:
            run(command, root / 'commands.jsonl', cwd=root, env=env,
                stdout=handle, stderr=subprocess.STDOUT)
        except BaseException:
            if (root / '.nextflow.log').is_file():
                save_nextflow_session(root)
            raise
    save_nextflow_session(root)
    records = rows(trace)
    expected = 2 if params['diagnostics'] else 1
    if len(records) != expected or any(r['status'] not in ('COMPLETED', 'CACHED') for r in records):
        raise ValueError(f'Native CelFiE execution did not complete: {records}')
    if require_cache and any(r['status'] != 'CACHED' for r in records):
        raise ValueError('Cached recovery unexpectedly recomputed native CelFiE')
    # Save the real native command, stdout and stderr alongside the final results.
    destination = root / 'native-task-evidence'
    destination.mkdir(exist_ok=True)
    import shutil
    for record in records:
        label = 'diagnostics' if 'OBSERVE_RESTARTS' in record['name'] else 'fit'
        workdir = Path(record['workdir'])
        task_fractions = workdir / '1_tissue_proportions.txt'
        published = root / ('diagnostics' if label == 'diagnostics' else 'results') / task_fractions.name
        if sha(task_fractions) != sha(published):
            raise ValueError('Published fractions differ from the completed native task output')
        if str(runtime / 'images/python.sif') not in (workdir / '.command.run').read_text():
            raise ValueError('Native task did not use the frozen local Python image')
        for name in ('.command.sh', '.command.out', '.command.err', '.command.run'):
            source = workdir / name
            target = destination / f'{label}{name}'
            if not target.exists():
                shutil.copy2(source, target)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--reference-counts', type=Path, required=True)
    parser.add_argument('--patient-counts', type=Path, required=True)
    parser.add_argument('--reference-validation', type=Path, required=True)
    parser.add_argument('--samples', nargs='+', required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--verify-cache', action='store_true')
    parser.add_argument('--diagnostics', action='store_true', help='Observe the same ten native starts in a separate replay')
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Submit run_celfie_patients.sh through sbatch')
    root = args.output.absolute()
    if not args.resume:
        root = new_run(root)
    else:
        require_incomplete(root)
    with run_lock(root):
        require_incomplete(root)
        preflights = root / 'preflight-attempts.jsonl'
        old_attempts = [json.loads(line)['attempt'] for line in preflights.read_text().splitlines()] if preflights.exists() else []
        record = dict(attempt=max(old_attempts, default=0) + 1, started_utc=utc(),
                      status='running', resume=args.resume, slurm_job_id=os.environ['SLURM_JOB_ID'])
        with preflights.open('a') as handle:
            handle.write(json.dumps(record) + '\n')
        try:
            prepared = prepare_run(args, root)
        except BaseException as exc:
            finish(root, 'failed', stage='preflight', error=str(exc))
            record.update(status='failed', ended_utc=utc(), error=str(exc))
            raise
        else:
            record.update(status='complete', completed_utc=utc())
        finally:
            with preflights.open('a') as handle:
                handle.write(json.dumps(record) + '\n')
        execute_run(args, root, *prepared)


def prepare_run(args, root):
    for name in ('nextflow', 'apptainer', 'sbatch', 'scontrol'):
        import shutil
        if not shutil.which(name):
            raise ValueError(f'Required command missing: {name}')
    runtime = verify_runtime(args.runtime)
    if not WORK.is_dir() or not os.access(WORK, os.W_OK):
        raise ValueError(f'Shared Nextflow work directory is missing or unwritable: {WORK}')
    reference_path = args.reference_counts.resolve(strict=True)
    patient_path = args.patient_counts.resolve(strict=True)
    gate_path = args.reference_validation.resolve(strict=True)
    atlas, cells = atlas_rows()
    coordinates = [(r['chr'], r['start'], r['end']) for r in atlas]
    reference, _ = count_table(reference_path, coordinates, cells)
    patients, available = count_table(patient_path, coordinates)
    if any(s not in available for s in args.samples):
        raise ValueError('Requested sample missing from patient count table')
    gate_reference(gate_path, reference_path, patient_path)
    merged = native_table(reference, patients, cells, args.samples)
    config = render_config(runtime, root, WORKFLOW / 'tools/celfie_patients.nf')
    params = dict(SETTINGS, deconflow_source=str(runtime / ('DecoNFlow-' + COMMIT)),
                  count_table=str(root / 'counts.tsv'), outdir=str(root / 'results'),
                  nsamples=len(args.samples), diagnostics=args.diagnostics,
                  diagnostic_script=str(WORKFLOW / 'tools/celfie_diagnostics.py'))
    run_identity = dict(schema_version=1, method='celfie', atlas=identity(ATLAS, digest=True),
                        reference_counts=identity(reference_path, digest=True),
                        patient_counts=identity(patient_path, digest=True),
                        reference_validation=identity(gate_path, digest=True),
                        runtime=str(runtime), runtime_lock=identity(runtime / 'runtime.lock.json', digest=True),
                        source_commit=COMMIT, settings=SETTINGS, initial_seed=42,
                        samples=args.samples, cells=cells, diagnostics=args.diagnostics,
                        code={p: sha(WORKFLOW / p) for p in OWN_FILES},
                        merged_counts_sha256=hashlib.sha256(merged.encode()).hexdigest(),
                        config_sha256=hashlib.sha256(config.encode()).hexdigest(),
                        versions=runtime_versions())
    return run_identity, runtime, params, merged, config, cells


def execute_run(args, root, run_identity, runtime, params, merged, config, cells):
    if args.resume:
        require_incomplete(root)
        if read_json(root / 'identity.json') != run_identity:
            raise ValueError('Refusing changed inputs, settings, source, runtime or configuration on resume')
        if not (root / 'execution.json').exists():
            raise ValueError('Run has no original execution-start timestamp')
    else:
        write_json(root / 'identity.json', run_identity)
    attempts = root / 'attempts.jsonl'
    previous_attempts = [json.loads(line)['attempt'] for line in attempts.read_text().splitlines()] if attempts.exists() else []
    attempt = max(previous_attempts, default=0) + 1
    attempt_record = dict(attempt=attempt, started_utc=utc(), status='running',
                          resume=args.resume, verify_cache=args.verify_cache,
                          slurm_job_id=os.environ['SLURM_JOB_ID'])
    with attempts.open('a') as handle:
        handle.write(json.dumps(attempt_record) + '\n')
    try:
        immutable_text(root / 'counts.tsv', merged)
        immutable_text(root / 'celfie.config', config)
        immutable_text(root / 'params.json', json.dumps(params, indent=2, allow_nan=False) + '\n')
        if {p: sha(WORKFLOW / p) for p in OWN_FILES} != run_identity['code']:
            raise ValueError('Local fitting code changed during run preparation')
        state = read_json(root / 'execution.json')
        state.update(status='running', methods=['celfie'], samples=args.samples,
                     upstream_unchanged=True, reference_validation=run_identity['reference_validation']['path'])
        write_json(root / 'execution.json', state)
        before = sha(root / 'results/1_tissue_proportions.txt') if (root / 'results/1_tissue_proportions.txt').exists() else None
        invoke(root, runtime, params, 'fit', attempt, resume=args.resume)
        result = validate_results(root, args.samples, cells, args.diagnostics)
        if before and result['fractions']['sha256'] != before:
            raise ValueError('Resume changed completed native fractions')
        if args.verify_cache:
            invoke(root, runtime, params, 'recovery', attempt, resume=True, require_cache=True)
            recovered = validate_results(root, args.samples, cells, args.diagnostics)
            if recovered['fractions']['sha256'] != result['fractions']['sha256']:
                raise ValueError('Cached recovery changed native fractions')
            result['cached_recovery_verified'] = True
        if {p: sha(WORKFLOW / p) for p in OWN_FILES} != run_identity['code']:
            raise ValueError('Local fitting code changed during native execution')
        complete = dict(result, completed_utc=utc(),
                        scope='Patient count integrity and native fraction output validation; restart convergence only when diagnostics requested',
                        biological_accuracy_validated=False,
                        missing_count_handling='Actual zero depth retained; native CelFiE pseudocount handling unchanged')
        write_json(root / 'COMPLETE.json', complete)
        write_json(root / 'VALIDATED.json', complete)
        finish(root)
        attempt_record.update(status='complete', completed_utc=utc())
    except Exception as exc:
        finish(root, 'failed', error=str(exc))
        attempt_record.update(status='failed', ended_utc=utc(), error=str(exc))
        raise
    finally:
        with attempts.open('a') as handle:
            handle.write(json.dumps(attempt_record) + '\n')


if __name__ == '__main__':
    main()
