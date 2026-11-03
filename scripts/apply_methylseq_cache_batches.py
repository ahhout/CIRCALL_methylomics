#!/usr/bin/env python3
"""Apply exact approved cache manifests together; dry-run unless --apply is set."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import io
import json
import os
from pathlib import Path
import stat
import tarfile
from unittest.mock import patch

import cleanup_methylseq_cache as c


def collect_parallel():
    scan = c.scan
    roots = [(r, True) for r in c.ROOTS] + [(c.REPO, False), (c.RESULTS, False),
             (c.BASE / 'downstream', False), (c.RESULTS.parent / 'references', False)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        pending = {(str(r), mode): pool.submit(scan, r, mode) for r, mode in roots}
        def cached(root, mode=False):
            key = (str(root), mode)
            return pending[key].result() if key in pending else scan(root, mode)
        with patch.object(c, 'scan', cached):
            return c.collect()


def load_manifests(paths, hashes):
    if not paths or len(paths) != len(hashes):
        raise c.Refused('Each manifest requires its exact approved hash')
    manifests, tasks, entries, batches = [], set(), set(), set()
    for path, digest in zip(paths, hashes):
        raw = Path(path).read_bytes()
        if c.sha(raw) != digest:
            raise c.Refused('Manifest approval hash mismatch')
        m = json.loads(raw)
        if m['blockers'] or m['helper_sha256'] != c.sha(Path(c.__file__).read_bytes()):
            raise c.Refused('Blocked manifest or changed original helper')
        if m['batch'] not in {'abandoned_alignment', 'superseded_integrated', 'legacy_emseq'} or c.atlas_concurrent(m['review']):
            raise c.Refused('Only idle integrated/legacy batches are supported')
        if m['batch'] in batches or not m['tasks'] or len(set(m['tasks'])) != len(m['tasks']):
            raise c.Refused('Duplicate or empty batch/task scope')
        for task in m['tasks']:
            if c.task_for(task) != task or task in tasks or not any(task.startswith(str(r) + '/') for r in c.ROOTS[:3]):
                raise c.Refused('Invalid or overlapping task scope')
            tasks.add(task)
        for row in m['entries']:
            if c.task_for(row['path']) not in m['tasks'] or row['path'] in entries:
                raise c.Refused('Invalid or overlapping entry scope')
            entries.add(row['path'])
        batches.add(m['batch'])
        manifests.append((str(path), digest, m))
    return manifests


def check_entries(m, state):
    fields = ('path', 'device', 'inode', 'mode', 'nlink', 'bytes', 'mtime_ns', 'link')
    def basic(rows):
        return sorted(({k: r[k] for k in fields if k in r} for r in rows), key=lambda r: r['path'])
    actual = [r for t in m['tasks'] for r in state['tasks'].get(t, [])]
    if basic(actual) != basic(m['entries']):
        raise c.Refused('Approved task contents changed')
    for task in m['tasks']:
        if Path(task).resolve(strict=True) != Path(task):
            raise c.Refused('Task path has symlink parents')


def preserve(entries, output):
    records = {}
    archive = output / 'task_provenance.tar'
    with archive.open('xb') as handle:
        with tarfile.open(fileobj=handle, mode='w') as bundle:
            for row in entries:
                p = Path(row['path'])
                if not (stat.S_ISREG(row['mode']) and row['bytes'] <= c.SMALL and p.name.startswith(('.command', '.exitcode'))):
                    continue
                data = p.read_bytes()
                current = c.identity(p)
                if any(current[k] != row[k] for k in current):
                    raise c.Refused('Provenance source changed')
                name = str(p.relative_to(c.BASE))
                member = tarfile.TarInfo(name); member.size = len(data); member.mode = 0o600
                bundle.addfile(member, io.BytesIO(data))
                records[name] = c.sha(data)
        handle.flush(); os.fsync(handle.fileno())
    with tarfile.open(archive) as bundle:
        members = bundle.getmembers()
        if len(members) != len(records) or {x.name for x in members} != set(records):
            raise c.Refused('Provenance archive membership mismatch')
        for member in members:
            if not member.isfile() or c.sha(bundle.extractfile(member).read()) != records[member.name]:
                raise c.Refused('Provenance archive verification failed')
    c.write_new(output / 'provenance_sha256.json', records)
    return len(records)


def delete_task(rows, handle):
    def log(event, **values):
        handle.write(json.dumps(dict(time=c.now(), event=event, **values)) + '\n')
        handle.flush()
    ordered = sorted(rows, key=lambda r: len(Path(r['path']).parts), reverse=True)
    # Durable intents precede every unlink; result records are synced per task.
    # A crash can leave an intent without a result, requiring reconciliation.
    for row in ordered:
        log('intent', path=row['path'])
    os.fsync(handle.fileno())
    try:
        for row in ordered:
            c.remove_entry(row)
            log('removed', path=row['path'])
    except BaseException as exc:
        log('stopped', reason=str(exc))
        raise
    finally:
        os.fsync(handle.fileno())


def run(paths, hashes, output, apply=False):
    manifests = load_manifests(paths, hashes)
    output = c.metadata_destination(output)
    output.mkdir(mode=0o700, parents=True)
    c.write_new(output / 'operation.json', {'started_utc': c.now(), 'apply': apply,
        'manifests': dict(zip(paths, hashes)), 'executor_sha256': c.sha(Path(__file__).read_bytes())})
    with c.maintenance():
        queue = c.scheduler()
        for _, _, m in manifests:
            c.require_idle(m['review'], queue)
        state = collect_parallel()
        for _, _, m in manifests:
            if c.unresolved_errors(state, m['review']) or set(c.batch_selection(state, m['batch'])) != set(m['tasks']):
                raise c.Refused('Dependency coverage or candidate scope changed')
            c.check_protected(m['state'], state, m['review'])
            if state['protected'] != m['state']['protected'] or state['findings'] != m['state']['findings'] or state['traces'] != m['state']['traces']:
                raise c.Refused('Dependency evidence changed')
            errors = c.validate_review(m['review'], state, m['tasks'])
            if errors:
                raise c.Refused('; '.join(errors))
            check_entries(m, state)
        c.write_new(output / 'preflight.json', {'checked_utc': c.now(), 'status': 'PASS', 'protected_tasks': len(state['protected'])})
        if not apply:
            return
        entries = [r for _, _, m in manifests for r in m['entries']]
        tasks = {t for _, _, m in manifests for t in m['tasks']}
        before = c.allocated_usage()
        c.write_new(output / 'allocated_before.json', before)
        count = preserve(entries, output)
        print(f'Preflight passed; {count} provenance files archived and verified.', flush=True)
        with (output / 'deletions.jsonl').open('x') as journal:
            for path, digest, m in manifests:
                c.require_idle(m['review'], c.scheduler())
                journal.write(json.dumps({'event': 'batch_start', 'time': c.now(), 'manifest': path, 'sha256': digest}) + '\n')
                journal.flush(); os.fsync(journal.fileno())
                grouped = {t: [] for t in m['tasks']}
                for row in m['entries']:
                    grouped[c.task_for(row['path'])].append(row)
                for i, (task, rows) in enumerate(grouped.items(), 1):
                    delete_task(rows, journal)
                    if i % 20 == 0 or i == len(grouped):
                        print(f'{m["batch"]}: {i}/{len(grouped)} tasks removed', flush=True)
            journal.write(json.dumps({'event': 'complete', 'time': c.now(), 'tasks': len(tasks)}) + '\n')
            journal.flush(); os.fsync(journal.fileno())
        if any(Path(t).exists() or Path(t).is_symlink() for t in tasks):
            raise c.Refused('Retired task remains')
        post = collect_parallel()
        for _, _, m in manifests:
            if c.unresolved_errors(post, m['review']):
                raise c.Refused('Post-delete dependency scan failed')
            c.check_protected(state, post, m['review'])
        links = lambda s: {(f['path'], f.get('target')) for f in s['findings'] if f['kind'] == 'external_link' and c.task_for(f['path']) not in tasks}
        if links(state) != links(post):
            raise c.Refused('Retained consumer links changed')
        events = [json.loads(line) for line in (output / 'deletions.jsonl').read_text().splitlines()]
        removed = [e['path'] for e in events if e['event'] == 'removed']
        if len(removed) != len(set(removed)) or set(removed) != {r['path'] for r in entries}:
            raise c.Refused('Deletion journal does not match exact approved entries')
        after = c.allocated_usage()
        c.write_new(output / 'execution.json', {'completed_utc': c.now(), 'status': 'PASS',
            'tasks': len(tasks), 'removed_entries': len(removed), 'protected_tasks_verified': len(state['protected']),
            'provenance_files_verified': count, 'allocated_before': before, 'allocated_after': after,
            'measured_reduction_bytes': sum(before.values()) - sum(after.values()),
            'manifest_hashes': dict(zip(paths, hashes))})
        print('Deletion and protected-file verification completed.', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--approved-manifest', action='append', required=True)
    p.add_argument('--confirm-manifest-sha256', action='append', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--apply', action='store_true')
    a = p.parse_args()
    run(a.approved_manifest, a.confirm_manifest_sha256, a.output_dir, a.apply)


if __name__ == '__main__':
    main()
