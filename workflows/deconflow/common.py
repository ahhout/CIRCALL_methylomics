"""Small shared provenance helpers; no scientific defaults live here."""
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from output_contract import output_guard

WORKFLOW = Path(__file__).resolve().parent
PROJECT = WORKFLOW.parents[1]
WORK = Path('/path/to/hpc-scratch/tmp/em_seq/work_emseq_integrated')
COMMIT = '4195f65ad986ea5a68129d2da5b2b1effed25bc8'
ATLAS = PROJECT / 'workflows/uxm/reference/Atlas.U25.l4.hg38.full.tsv'
ATLAS_SHA = 'e58fbec4fe6200a311540c6e1a16be901a2d1e9ed9d74b6a57888f2feb635f51'
GENOME = (PROJECT / 'references/wgbstools_hg38').resolve()
CPG_SHA = '28180d18667dc0a6e66dd4250593e03a8bf10be50e9212aefd4e9282afea8590'


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def rows(path, delimiter='\t'):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def table(path, records, columns, delimiter='\t'):
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, columns, delimiter=delimiter, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


def new_run(path):
    path = Path(path).absolute()
    if not re.fullmatch(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-.+-v[1-9]\d*', path.name):
        raise ValueError('New bundle needs YYYY-MM-DD_HH-MM-<purpose>-vN in UTC')
    dt.datetime.strptime(path.name[:16], '%Y-%m-%d_%H-%M')
    output_guard.reserve_path(path, {'workflow': str(WORKFLOW), 'code_sha256': sha(__file__),
                                     'project_revision': output_guard.require_clean_source(PROJECT)})
    output_guard.update(path, status='running')
    write_json(path / 'execution.json', {'started_utc': utc(), 'status': 'running'})
    snapshot = path / 'provenance/code'
    snapshot.mkdir(parents=True)
    if os.environ.get('SLURM_JOB_ID'):
        write_json(path / 'provenance/slurm_environment.json',
                   {k: v for k, v in os.environ.items() if k.startswith('SLURM_')})
        allocation = subprocess.run(['scontrol', 'show', 'job', os.environ['SLURM_JOB_ID']],
                                    capture_output=True, text=True,
                                    env=dict(os.environ, TZ='UTC'), check=True)
        (path / 'provenance/slurm_job_utc.txt').write_text(allocation.stdout)
    for directory in (WORKFLOW, WORKFLOW / 'tools', WORKFLOW / 'config'):
        for code in sorted(p for p in directory.iterdir() if p.suffix in ('.py', '.sh', '.nf', '.config') or p.name.startswith('params.')):
            target = snapshot / code.relative_to(WORKFLOW)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(code, target)
    write_json(path / 'provenance/code.sha256.json', {str(p.relative_to(snapshot)): sha(p) for p in snapshot.rglob('*') if p.is_file()})
    return path


def finish(path, status='complete', **fields):
    p = Path(path) / 'execution.json'
    record = read_json(p)
    record.update(status=status, **fields)
    record['completed_utc' if status == 'complete' else 'ended_utc'] = utc()
    write_json(p, record)
    output_guard.finish(path, status == 'complete')


def require_incomplete(path):
    revision = output_guard.require_clean_source(PROJECT)
    record = output_guard.read_record(path)
    if record['status'] not in ('running', 'failed'):
        raise ValueError('Only owned incomplete runs may resume')
    if record['identity'].get('project_revision') != revision:
        raise ValueError('Resume requires the original committed project revision')
    if record['identity']['code_sha256'] != sha(__file__):
        raise ValueError('Resume requires original shared helper source')
    for marker in ('COHORT_COMPLETE', 'COMPLETE.json', 'COMPLETE', 'VALIDATED.json'):
        if (Path(path) / marker).exists():
            raise ValueError('Completed outputs are immutable')


def run(command, log=None, **kwargs):
    if log:
        with Path(log).open('a') as handle:
            handle.write(json.dumps({'started_utc': utc(), 'argv': [str(x) for x in command]}) + '\n')
    return subprocess.run([str(x) for x in command], check=True, **kwargs)


def identity(path, digest=False):
    p = Path(path).resolve(strict=True)
    s = p.stat()
    result = {'path': str(p), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}
    if digest:
        result['sha256'] = sha(p)
    return result


def atlas_rows():
    if sha(ATLAS) != ATLAS_SHA:
        raise ValueError('Canonical UXM atlas checksum changed')
    data = rows(ATLAS)
    cells = list(data[0])[8:]
    if len(cells) != 40 or len({r['name'] for r in data}) != len(data):
        raise ValueError('Unexpected atlas columns or duplicate marker names')
    return data, cells


def current_nextflow_session(launch):
    """Resolve exactly the logged session and verify its local cache/history."""
    launch = Path(launch)
    log = launch / '.nextflow.log'
    if log.is_symlink() or not log.is_file():
        raise ValueError('Missing local Nextflow session log; cannot resume')
    matches = set(re.findall(r'Session UUID:\s*([0-9a-f]{8}-(?:[0-9a-f]{4}-){3}[0-9a-f]{12})', log.read_text()))
    if len(matches) != 1:
        raise ValueError('Missing or conflicting Nextflow session UUIDs')
    session = matches.pop()
    cache = launch / '.nextflow/cache' / session / 'db'
    history = launch / '.nextflow/history'
    if cache.is_symlink() or not cache.is_dir() or not history.is_file() or session not in history.read_text():
        raise ValueError('Nextflow session cache/history missing or inconsistent')
    return session


def save_nextflow_session(launch):
    session = current_nextflow_session(launch)
    receipt = Path(launch) / '.nextflow-session.json'
    record = {'session_id': session, 'launch_directory': str(Path(launch).resolve())}
    if receipt.exists() and read_json(receipt) != record:
        raise ValueError('Nextflow session changed within an existing execution')
    if not receipt.exists():
        with receipt.open('x') as handle:
            json.dump(record, handle, indent=2)
    return session


def resume_nextflow_session(launch):
    receipt = Path(launch) / '.nextflow-session.json'
    if receipt.is_symlink() or not receipt.is_file():
        raise ValueError('Missing original Nextflow session receipt')
    record = read_json(receipt)
    if record != {'session_id': current_nextflow_session(launch), 'launch_directory': str(Path(launch).resolve())}:
        raise ValueError('Nextflow session identity differs from original execution')
    return record['session_id']
