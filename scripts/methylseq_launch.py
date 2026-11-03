#!/usr/bin/env python3
"""Launch pinned methylseq with isolated provenance and exclusive output ownership.

The original methylseq_cache_session.py remains a historical recovery interface.
This launcher never reads its legacy UUID or the original runtime database.
"""
import argparse
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import uuid

sys.path.insert(0, str(Path(__file__).resolve().parent / 'lib'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import output_guard

PROJECT = Path(__file__).resolve().parents[1]
NAME = re.compile(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-emseq-integrated-[a-z0-9]+(?:-[a-z0-9]+)*-v[1-9][0-9]*')
UUID = re.compile(r'[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}')
PATH_KEYS = {'pipeline_root', 'runtime_root', 'work_root', 'maintenance_lock',
             'results_parent', 'reference_root', 'container_cache'}


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(value, message):
    if not value:
        raise ValueError(message)


def git(repo, *args):
    result = subprocess.run(['git', '--no-optional-locks', '-C', str(repo), *args],
                            check=True, capture_output=True)
    return result.stdout.decode().strip()


def clean_commit(repo):
    require(Path(git(repo, 'rev-parse', '--show-toplevel')) == repo,
            'Source must be the canonical Git checkout root')
    require(not git(repo, 'status', '--porcelain', '--untracked-files=all'),
            f'Commit reviewed source changes before execution: {repo}')
    return git(repo, 'rev-parse', 'HEAD')


def load_paths(project):
    raw = json.loads((project / 'config/methylseq_paths.json').read_text())
    require(set(raw) == PATH_KEYS | {'schema'} and raw['schema'] == 1,
            'Invalid methylseq path configuration')
    paths = {key: Path(raw[key]) for key in PATH_KEYS}
    for key, path in paths.items():
        require(path.is_absolute() and '..' not in path.parts, f'Invalid absolute path: {key}')
        if key != 'reference_root':
            require(path == path.resolve(), f'Symlink or noncanonical path rejected: {key}')
    require(paths['runtime_root'] == project / 'runs/methylseq',
            "New runtime must belong to this project's runs/methylseq directory")
    for key in ('pipeline_root', 'work_root', 'results_parent', 'reference_root', 'container_cache'):
        require(paths[key].is_dir(), f'Missing existing configured directory: {key}')
    require(paths['maintenance_lock'].parent.is_dir(), 'Missing shared lock parent')
    return paths


def input_files(project, paths, smoke_inputs=None):
    files = {
        'run_methylseq_integrated.sh': project / 'run_methylseq_integrated.sh',
        'params.emseq.integrated.yml': project / 'config/params.emseq.integrated.yml',
        'nextflow.config': project / 'config/nextflow.config',
        'samplesheet.csv': project / 'samplesheets/samplesheet.csv',
        'reference_manifest.tsv': paths['reference_root'] / 'manifest.tsv',
        'methylseq_paths.json': project / 'config/methylseq_paths.json',
        'methylseq_pipeline.json': project / 'config/methylseq_pipeline.json',
    }
    if smoke_inputs is not None:
        import methylseq_smoke
        manifest = methylseq_smoke.validate(project, Path(smoke_inputs))
        files['production_samplesheet.csv'] = files['samplesheet.csv']
        files['samplesheet.csv'] = Path(manifest['samplesheet'])
        files['smoke_inputs.json'] = Path(smoke_inputs)
        files['smoke_resources.config'] = project / 'config/methylseq_smoke_resources.config'
    return files


def flat_params(path):
    """Read the project's deliberately flat scalar YAML without a new dependency."""
    values = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        match = re.fullmatch(r'([a-zA-Z_][a-zA-Z_0-9]*):\s*(.*?)\s*', line)
        require(match is not None, 'Production params must remain a flat scalar YAML mapping')
        key, value = match.groups()
        require(key not in values, 'Duplicate production parameter')
        require(value and not value.startswith(('[', '{', '&', '*', '!', '|', '>')),
                'Complex YAML values require an explicit launcher update')
        values[key] = value.strip('\"\'')
    return values


def check_publication(project, paths):
    params = flat_params(project / 'config/params.emseq.integrated.yml')
    forbidden = {'outdir', 'custom_config_base', 'custom_config_version', 'trace_report_suffix'}
    require(not forbidden.intersection(params), 'Conflicting output or external configuration override in params')
    expected = {'input': project / 'samplesheets/samplesheet.csv',
                'fasta': paths['reference_root'] / 'BismarkIndex/genome.fa',
                'bismark_index': paths['reference_root'] / 'BismarkIndex'}
    for key, path in expected.items():
        require(params.get(key) == str(path), f'Configured {key} differs from project path contract')
    config = (project / 'config/nextflow.config').read_text()
    uncommented = re.sub(r'/\*.*?\*/|//[^\n]*', '', config, flags=re.S)
    require(not re.search(r'\b(?:outdir|publishDir|includeConfig|timeline|report|trace|dag)\b', uncommented),
            'Project config must not override publication or include external configuration')
    require(not any(os.environ.get(key) for key in ('NXF_CONFIG', 'NXF_PARAMS', 'NXF_CLI')),
            'Unset external Nextflow configuration/parameter overrides')
    require(not (paths['runtime_root'] / 'nextflow.config').exists(),
            'Unexpected runtime nextflow.config; use the committed project config')


def preflight(project, run_name, smoke_inputs=None):
    project = Path(project).resolve(strict=True)
    require(NAME.fullmatch(run_name), 'Use YYYY-MM-DD_HH-MM-emseq-integrated-<purpose>-vN (UTC)')
    datetime.strptime(run_name[:16], '%Y-%m-%d_%H-%M')
    paths = load_paths(project)
    require((project / 'scripts/summarize_emseq_controls.py').is_file(), 'Missing control-QC reporter')
    files = input_files(project, paths, smoke_inputs)
    for name, path in files.items():
        require(path.is_file() and path.stat().st_size, f'Missing or empty required input: {name}')
    check_publication(project, paths)
    reference = paths['reference_root'] / 'BismarkIndex'
    require((reference / 'genome.fa').is_file(), 'Missing combined reference FASTA')
    for conversion in ('CT_conversion', 'GA_conversion'):
        directory = reference / 'Bisulfite_Genome' / conversion
        require(directory.is_dir(), 'Missing Bismark reference conversion directory')
        require(any(directory.glob('*.bt2')) or any(directory.glob('*.bt2l')), 'Missing Bowtie2 index files')
    # Only presence checks: no identifiers in console output or bulk data hashing.
    with files['samplesheet.csv'].open(newline='') as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ['sample', 'fastq_1', 'fastq_2', 'genome'], 'Unexpected production samplesheet columns')
        count = 0
        for row in reader:
            require(None not in row and row['sample'] and row['fastq_1'], 'Malformed samplesheet row')
            for field in ('fastq_1', 'fastq_2'):
                if row[field]:
                    require(Path(row[field]).is_file() and Path(row[field]).stat().st_size,
                            'A samplesheet FASTQ is missing or empty')
            count += 1
        require(count, 'Empty samplesheet')
    pipeline = paths['pipeline_root']
    pin = json.loads(files['methylseq_pipeline.json'].read_text())
    require(pin.get('schema') == 1 and pin.get('pipeline') == 'nf-core/methylseq'
            and pin.get('release') == '4.2.0', 'Invalid pipeline lock')
    require((pipeline / 'main.nf').is_file() and (pipeline / 'nextflow.config').is_file(), 'Pipeline source is incomplete')
    pipeline_commit = clean_commit(pipeline)
    require(pipeline_commit == pin['commit'], 'Pipeline checkout differs from pinned commit')
    for key in ('upstream_release_commit', 'fork_baseline_commit'):
        git(pipeline, 'merge-base', '--is-ancestor', pin[key], pipeline_commit)
    identity = {'schema': 1, 'project_source': str(project), 'project_commit': clean_commit(project),
                'pipeline_source': str(pipeline), 'pipeline_commit': pipeline_commit,
                'pipeline_release': '4.2.0', 'runtime_root': str(paths['runtime_root']),
                'work_root': str(paths['work_root']),
                'inputs': {name: sha(path) for name, path in files.items()}}
    return paths, files, identity


@contextmanager
def execution_lock(path):
    # Same inode as existing canonical launchers and cache maintenance.
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | getattr(os, 'O_NOFOLLOW', 0), 0o600)
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError('Shared cache is in use or under maintenance') from exc
        yield
    finally:
        os.close(descriptor)


def option(words, name):
    if words.count(name) != 1:
        return None
    index = words.index(name)
    return words[index + 1] if index + 1 < len(words) else None


def output_session(runtime, output, paths):
    history = runtime / '.nextflow/history'
    require(history.is_file(), 'No new-layout Nextflow session history is available')
    sessions = set()
    for line in history.read_text().splitlines():
        fields = line.split('\t', 6)
        if len(fields) != 7:
            continue
        words = shlex.split(fields[6])
        if option(words, '--outdir') != str(output):
            continue
        require(option(words, '-work-dir') == str(paths['work_root']), 'Session work-root mismatch')
        require(option(words, 'run') == str(paths['pipeline_root']), 'Session pipeline-source mismatch')
        require(UUID.fullmatch(fields[5]), 'Invalid recorded session UUID')
        sessions.add(fields[5])
    require(len(sessions) == 1, 'Missing or conflicting session mapping for this exact output')
    session = sessions.pop()
    require((runtime / '.nextflow/cache' / session / 'db/CURRENT').is_file(), 'Recorded session cache metadata is missing')
    return session


def verify_snapshots(submitted, files, identity):
    require(not submitted.is_symlink() and not submitted.parent.is_symlink(),
            'Submitted snapshot directories must not be symlinks')
    recorded = json.loads((submitted / 'identity.json').read_text())
    require(recorded == identity, 'Submitted source/input identity changed')
    for name, path in files.items():
        snapshot = submitted / name
        require(not snapshot.is_symlink() and snapshot.is_file()
                and sha(snapshot) == identity['inputs'][name] == sha(path),
                f'Submitted input snapshot changed: {name}')
    expected_checks = ''.join(f'{digest}  {name}\n' for name, digest in sorted(identity['inputs'].items()))
    require((submitted / 'SHA256SUMS').read_text() == expected_checks,
            'Submitted checksum inventory changed')
    record = output_guard.read_record(submitted.parent.parent)
    require(sha(submitted / 'run_manifest.json') == record.get('submitted_manifest_sha256'),
            'Original run manifest changed')


def check_destination(paths, files, identity, run_name, resume_existing):
    """Read-only advice; execution repeats ownership checks while holding the lock."""
    output = paths['results_parent'] / run_name
    if resume_existing:
        record = output_guard.resume(output, identity)
        session = record.get('session_id')
        require(session and UUID.fullmatch(session), 'Incomplete run lacks a recorded session; do not adopt it')
        require(output_session(paths['runtime_root'], output, paths) == session,
                'Run session mapping differs from history')
        verify_snapshots(output / 'pipeline_info/submitted', files, identity)
    else:
        require(not os.path.lexists(output), 'Output already exists; choose a new run name')
    link = Path(identity['project_source']) / 'outputs' / run_name
    if os.path.lexists(link):
        require(link.is_symlink() and link.resolve() == output,
                'Navigation path already belongs to another output')
    return output


def exclusive_json(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write('\n')


def groovy_string(path):
    return "'" + str(path).replace('\\', '\\\\').replace("'", "\\'") + "'"


def command_for(paths, output, submitted, attempt, session=None):
    configs = [paths['pipeline_root'] / 'nextflow.config', submitted / 'nextflow.config', attempt / 'publication.config']
    if (submitted / 'smoke_resources.config').is_file():
        configs.insert(2, submitted / 'smoke_resources.config')
    command = ['nextflow', '-C', ','.join(map(str, configs)), '-log', str(attempt / 'nextflow.log'),
               'run', str(paths['pipeline_root']), '-params-file', str(submitted / 'params.emseq.integrated.yml'),
               '-work-dir', str(paths['work_root']), '-profile', 'apptainer',
               '--custom_config_base', 'false',
               '--input', str(submitted / 'samplesheet.csv'), '--outdir', str(output)]
    if session:
        command += ['-resume', session]
    return command


def check_versions():
    nxf = subprocess.run(['nextflow', '-version'], check=True, text=True, capture_output=True)
    require(re.search(r'\bversion 25\.10\.2\b', nxf.stdout + nxf.stderr), 'Expected Nextflow 25.10.2')
    apptainer = subprocess.check_output(['apptainer', '--version'], text=True).strip()
    require(apptainer == 'apptainer version 1.3.6', 'Expected module-provided Apptainer 1.3.6')
    return {'nextflow': (nxf.stdout + nxf.stderr).strip(), 'nextflow_executable': shutil.which('nextflow'),
            'apptainer': apptainer, 'apptainer_executable': shutil.which('apptainer'),
            'apptainer_module': 'containers/apptainer/1.3.6', 'conda_environment': 'nextflow'}


def execute(project, run_name, resume_existing=False, *, runner=subprocess.run, smoke_inputs=None):
    project = Path(project).resolve(strict=True)
    job = os.environ.get('SLURM_JOB_ID', '')
    require(re.fullmatch(r'[0-9]+', job), 'Submit a full execution with sbatch')
    paths, files, identity = preflight(project, run_name, smoke_inputs)
    versions = check_versions()
    output = paths['results_parent'] / run_name
    with execution_lock(paths['maintenance_lock']):
        require(preflight(project, run_name, smoke_inputs)[2] == identity, 'Source or input changed during preflight')
        runtime = paths['runtime_root']
        submitted = output / 'pipeline_info/submitted'
        session = None
        if resume_existing:
            record = output_guard.resume(output, identity)
            session = record.get('session_id')
            require(session and UUID.fullmatch(session), 'Incomplete run lacks a recorded session; do not adopt it')
            require(output_session(runtime, output, paths) == session, 'Run session mapping differs from history')
            verify_snapshots(submitted, files, identity)
        else:
            history = runtime / '.nextflow/history'
            if history.is_file():
                for line in history.read_text().splitlines():
                    fields = line.split('\t', 6)
                    if len(fields) == 7:
                        require(option(shlex.split(fields[6]), '--outdir') != str(output),
                                'Run name already occurs in runtime history; choose a new name')
            output = output_guard.reserve(paths['results_parent'], run_name, identity)
        attempt = None
        started = utc()
        try:
            if not resume_existing:
                submitted.mkdir(parents=True)
                for name, path in files.items():
                    with (submitted / name).open('xb') as handle:
                        handle.write(path.read_bytes())
                exclusive_json(submitted / 'identity.json', identity)
                with (submitted / 'SHA256SUMS').open('x') as handle:
                    handle.writelines(f'{digest}  {name}\n' for name, digest in sorted(identity['inputs'].items()))
                exclusive_json(submitted / 'run_manifest.json', {
                    'run_id': run_name, 'started_utc': started,
                    'run_named_utc': datetime.strptime(run_name[:16], '%Y-%m-%d_%H-%M').replace(tzinfo=timezone.utc).isoformat(),
                    'output': str(output), 'identity': identity, 'initial_resume': False})
                output_guard.update(output, submitted_manifest_sha256=sha(submitted / 'run_manifest.json'))
            verify_snapshots(submitted, files, identity)
            runtime.mkdir(parents=True, exist_ok=True)
            attempt_id = job + '-' + uuid.uuid4().hex[:12]
            attempt = submitted / 'attempts' / attempt_id
            attempt.mkdir(parents=True, exist_ok=False)
            with (attempt / 'publication.config').open('x') as handle:
                for kind, filename in (('timeline', 'timeline.html'), ('report', 'report.html'),
                                       ('trace', 'trace.txt'), ('dag', 'dag.html')):
                    handle.write(f'{kind}.file = {groovy_string(attempt / filename)}\n')
            command = command_for(paths, output, submitted, attempt, session)
            provenance = {'attempt_id': attempt_id, 'job_id': job, 'started_utc': started,
                          'status': 'running', 'session_id': session, 'resume_existing': resume_existing,
                          'command': command, 'identity': identity, 'runtime_environment': versions}
            exclusive_json(attempt / 'started.json', provenance)
            with (attempt / 'command.sh').open('x') as handle:
                handle.write(shlex.join(command) + '\n')
            shutil.copyfile(files['run_methylseq_integrated.sh'], attempt / 'launcher.sh')
            record = output_guard.read_record(output)
            output_guard.update(output, status='running', attempts=record.get('attempts', []) + [provenance])
            navigation = project / 'outputs'
            require(not navigation.is_symlink(), 'Navigation parent must not be a symlink')
            navigation.mkdir(exist_ok=True)
            output_guard.publish_link(navigation / run_name, output)
            environment = dict(os.environ, NXF_APPTAINER_CACHEDIR=str(paths['container_cache']),
                               APPTAINER_CACHEDIR=str(paths['container_cache']), CIRCALL_CONTAINER_CACHE=str(paths['container_cache']))
            require(preflight(project, run_name, smoke_inputs)[2] == identity, 'Source or input changed immediately before execution')
            result = runner(command, cwd=runtime, env=environment, check=False)
            try:
                observed = output_session(runtime, output, paths)
                require(session is None or session == observed, 'Pipeline changed session during resume')
                session = observed
                output_guard.update(output, session_id=session)
            except (OSError, ValueError):
                if result.returncode == 0:
                    raise
            require(result.returncode == 0, f'Nextflow exited with status {result.returncode}')
            control = output / 'control_qc/control_qc.tsv'
            require(not control.parent.is_symlink() and not control.is_symlink(),
                    'Control-QC publication must not follow symlinks')
            result = runner([sys.executable, str(project / 'scripts/summarize_emseq_controls.py'),
                             '--coverage-dir', str(output / 'bismark/methylation_calls/methylation_coverage'),
                             '--samplesheet', str(submitted / 'samplesheet.csv'), '--output', str(control),
                             '--owned-output-root', str(output), '--attempt-id', attempt_id],
                            cwd=runtime, env=environment, check=False)
            require(result.returncode == 0, f'Control QC exited with status {result.returncode}')
            require(control.is_file(), 'Control QC did not create the expected result')
            exclusive_json(attempt / 'ended.json', {'status': 'completed', 'completed_utc': utc(),
                                                   'exit_status': 0, 'session_id': session})
            output_guard.finish(output, True, session_id=session)
        except BaseException as exc:
            if session is None and attempt is not None:
                try:
                    session = output_session(runtime, output, paths)
                    output_guard.update(output, session_id=session)
                except (OSError, ValueError):
                    pass
            if attempt is not None and not (attempt / 'ended.json').exists():
                exclusive_json(attempt / 'ended.json', {'status': 'failed_or_interrupted', 'ended_utc': utc(),
                                                       'error': str(exc), 'session_id': session})
            output_guard.finish(output, False, session_id=session, last_error=str(exc))
            raise
    return output


def submission_command(project, run_name, resume_existing, identity, smoke_inputs=None):
    log_root = project / 'runs/methylseq/scheduler' / (run_name + '-' + uuid.uuid4().hex[:12])
    export = ','.join(['ALL', f'CIRCALL_PROJECT_ROOT={project}',
                       f"CIRCALL_PROJECT_COMMIT={identity['project_commit']}",
                       f"CIRCALL_PIPELINE_COMMIT={identity['pipeline_commit']}",
                       f'RUN_NAME={run_name}', f'RESUME_EXISTING={str(resume_existing).lower()}'])
    export += ',CIRCALL_SMOKE_INPUTS=' + (str(smoke_inputs) if smoke_inputs else '')
    resources = ['--partition=shortq', '--time=06:00:00'] if smoke_inputs else []
    return log_root, ['sbatch', '--export=' + export, '--open-mode=append',
                      '--output=' + str(log_root / '%j.out'), '--error=' + str(log_root / '%j.err'),
                      *resources, str(project / 'run_methylseq_integrated.sh')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['preflight', 'submit', 'run'])
    parser.add_argument('--project-root', type=Path, default=PROJECT)
    parser.add_argument('--run-name', default=os.environ.get('RUN_NAME'))
    parser.add_argument('--resume-existing', action='store_true', default=os.environ.get('RESUME_EXISTING', 'false') == 'true')
    parser.add_argument('--smoke-inputs', type=Path, default=os.environ.get('CIRCALL_SMOKE_INPUTS') or None,
                        help='Verified bounded paired-read manifest created by methylseq_smoke.py')
    args = parser.parse_args()
    require(os.environ.get('RESUME_EXISTING', 'false') in {'true', 'false'}, 'RESUME_EXISTING must be true or false')
    require(args.run_name, 'RUN_NAME or --run-name is required')
    project = args.project_root.resolve(strict=True)
    paths, files, identity = preflight(project, args.run_name, args.smoke_inputs)
    if args.action in {'preflight', 'submit'}:
        check_destination(paths, files, identity, args.run_name, args.resume_existing)
    if args.action == 'preflight':
        print(json.dumps({'status': 'preflight_passed', 'identity': identity,
                          'output': str(paths['results_parent'] / args.run_name),
                          'resume_existing': args.resume_existing,
                          'limits': 'No Nextflow execution, output reservation, or cache writes'}, indent=2))
    elif args.action == 'submit':
        log_root, command = submission_command(project, args.run_name, args.resume_existing, identity, args.smoke_inputs)
        log_root.mkdir(parents=True, exist_ok=False)
        subprocess.run(command, check=True)
    else:
        for variable, key in (('CIRCALL_PROJECT_COMMIT', 'project_commit'), ('CIRCALL_PIPELINE_COMMIT', 'pipeline_commit')):
            require(os.environ.get(variable) == identity[key], 'Submission commit missing or changed; use launcher --submit')
        print(execute(project, args.run_name, args.resume_existing, smoke_inputs=args.smoke_inputs))


if __name__ == '__main__':
    def interrupted(signum, frame):
        raise KeyboardInterrupt(f'Received signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, KeyboardInterrupt) as exc:
        sys.exit(f'Methylseq launch refused or incomplete: {exc}')
