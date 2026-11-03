#!/usr/bin/env python3
"""Validate explicit shared-session selection and record launch provenance."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import sys
import tempfile

REPO = Path(__file__).resolve().parents[1]
RECORD = REPO / 'config/methylseq_cache_session.json'
EXPECTED_WORK = '/path/to/hpc-scratch/tmp/em_seq/work_emseq_integrated'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write('\n')
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def validate_session(record=RECORD, repo=REPO):
    value = json.loads(Path(record).read_text())
    if set(value) != {'schema', 'session_id', 'work_dir'} or value['schema'] != 1:
        raise ValueError('Invalid shared session record')
    sid = value['session_id']
    if not re.fullmatch(r'[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}', sid):
        raise ValueError('Invalid session UUID')
    if value['work_dir'] != EXPECTED_WORK:
        raise ValueError('Shared work root differs from approved integrated root')
    if not (repo / '.nextflow/cache' / sid / 'db/CURRENT').is_file():
        raise ValueError('Recorded session cache metadata is missing')
    sessions = set()
    for line in (repo / '.nextflow/history').read_text().splitlines():
        fields = line.split('\t', 6)
        if len(fields) != 7:
            continue
        words = shlex.split(fields[6])
        if '-work-dir' in words and words[words.index('-work-dir') + 1] == value['work_dir']:
            sessions.add(fields[5])
    if sessions != {sid}:
        raise ValueError(f'Missing/conflicting integrated-session history: {sorted(sessions)}')
    return value


def verify_snapshots(submitted, repo=REPO):
    """Verify original snapshots, then compare scientific inputs; launcher may evolve."""
    names = {'run_methylseq_integrated.sh', 'params.emseq.integrated.yml',
             'nextflow.config', 'samplesheet.csv', 'reference_manifest.tsv'}
    checks = {}
    for line in (submitted / 'SHA256SUMS').read_text().splitlines():
        sha, name = line.split(None, 1)
        name = name.lstrip('*')
        if name not in names or name in checks or not re.fullmatch('[0-9a-f]{64}', sha):
            raise ValueError('Invalid submitted snapshot checksum inventory')
        checks[name] = sha
    if set(checks) != names:
        raise ValueError('Incomplete submitted snapshot checksum inventory')
    for name, sha in checks.items():
        if (submitted / name).is_symlink() or digest(submitted / name) != sha:
            raise ValueError(f'Original submitted snapshot changed: {name}')
    current = {'params.emseq.integrated.yml': repo / 'config/params.emseq.integrated.yml',
               'nextflow.config': repo / 'config/nextflow.config',
               'samplesheet.csv': repo / 'samplesheets/samplesheet.csv',
               'reference_manifest.tsv': repo / 'references/emseq_hg38_lambda_puc19/manifest.tsv'}
    for name, path in current.items():
        if digest(path) != checks[name]:
            raise ValueError(f'Current scientific input differs from submitted snapshot: {name}')


def infer_output_session(output, repo=REPO):
    ids = set()
    for line in (repo / '.nextflow/history').read_text().splitlines():
        parts = line.split('\t', 6)
        if len(parts) != 7:
            continue
        words = shlex.split(parts[6])
        if '--outdir' in words and words[words.index('--outdir') + 1] == str(output):
            ids.add(parts[5])
    if len(ids) != 1:
        raise ValueError('Existing run lacks one unambiguous historical session mapping')
    return ids.pop()


def begin(output, job, existing, started_utc=None):
    started_utc = started_utc or datetime.now(timezone.utc).isoformat()
    stamp = datetime.fromisoformat(started_utc.replace("Z", "+00:00"))
    if stamp.utcoffset() is None or stamp.utcoffset().total_seconds() != 0:
        raise ValueError("Attempt start must be timezone-labelled UTC")
    session = validate_session()
    submitted = output / 'pipeline_info/submitted'
    mapping = submitted / 'cache-session.json'
    if existing:
        verify_snapshots(submitted)
        recorded = json.loads(mapping.read_text()) if mapping.exists() else None
        sid = recorded['session_id'] if recorded else infer_output_session(output)
        if sid != session['session_id'] or (recorded and recorded != session):
            raise ValueError('Named run conflicts with shared session record')
    if mapping.exists() and json.loads(mapping.read_text()) != session:
        raise ValueError('Existing session mapping differs')
    if not mapping.exists():
        atomic_json(mapping, session)
    attempt = submitted / 'attempts' / job
    attempt.mkdir(parents=True, exist_ok=False)
    launcher = REPO / 'run_methylseq_integrated.sh'
    (attempt / 'launcher.sh').write_bytes(launcher.read_bytes())
    atomic_json(attempt / 'provenance.json', dict(session, output=str(output), job_id=job,
                launcher_sha256=digest(launcher), session_helper_sha256=digest(__file__),
                status='started', started_utc=started_utc, membership_complete=False))


def finish(output, job, status):
    attempt = output / 'pipeline_info/submitted/attempts' / job
    provenance = json.loads((attempt / 'provenance.json').read_text())
    traces = []
    tasks = []
    # Preserve available traces, including earlier attempts: absence is explicit.
    for trace in sorted((output / 'pipeline_info').glob('*trace*.txt')):
        data = trace.read_bytes()
        (attempt / trace.name).write_bytes(data)
        traces.append(dict(path=str(trace), sha256=hashlib.sha256(data).hexdigest()))
        with trace.open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                tasks.append({k: row.get(k) for k in ('hash', 'name', 'status', 'native_id')})
    ended_utc = datetime.now(timezone.utc).isoformat()
    provenance['ended_utc'] = ended_utc
    if status == 0:
        provenance['completed_utc'] = ended_utc
    else:
        provenance.pop('completed_utc', None)
    provenance.update(status='completed' if status == 0 else 'failed_or_interrupted',
                      exit_status=status, membership_complete=status == 0 and bool(traces), traces=traces)
    atomic_json(attempt / 'tasks.json', tasks)
    atomic_json(attempt / 'provenance.json', provenance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['validate', 'begin', 'finish'])
    parser.add_argument('--output', type=Path)
    parser.add_argument('--job')
    parser.add_argument('--existing', action='store_true')
    parser.add_argument('--exit-status', type=int)
    parser.add_argument('--started-utc')
    args = parser.parse_args()
    if args.action == 'validate':
        print(validate_session()['session_id'])
    else:
        if not args.output or not args.job or not re.fullmatch('[0-9]+', args.job):
            parser.error('output and numeric job are required')
        if args.action == 'begin':
            begin(args.output, args.job, args.existing, args.started_utc)
        else:
            if args.exit_status is None:
                parser.error('finish requires --exit-status')
            finish(args.output, args.job, args.exit_status)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        sys.exit(f'Session provenance error: {exc}')
