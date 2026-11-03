#!/usr/bin/env python3
"""Reconcile and resume an interrupted, approved atlas-cache deletion."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import stat
import subprocess
from unittest.mock import patch

import cleanup_methylseq_cache as cache


def mapped_devices(value, mounts, to_original=True):
    if isinstance(value, list):
        return [mapped_devices(v, mounts, to_original) for v in value]
    if not isinstance(value, dict):
        return value
    result = {k: mapped_devices(v, mounts, to_original) for k, v in value.items()}
    if 'path' in result and 'device' in result:
        for prefix, old, current in mounts:
            if result['path'].startswith(prefix + '/'):
                source, target = (current, old) if to_original else (old, current)
                if result['device'] != source:
                    raise cache.Refused('Unexpected filesystem device: ' + result['path'])
                result['device'] = target
                break
    return result


def mount_mapping(manifest):
    mapping = []
    rows = manifest['entries'] + manifest['state']['findings']
    for prefix, source, filesystem in (('/path/to/hpc-storage', 'beegfs_nodev', 'beegfs'),
                                       ('/home', 'nfs02:/zfs/home', 'nfs4')):
        result = subprocess.run(['findmnt', '-T', prefix, '-n', '-o', 'TARGET,SOURCE,FSTYPE'],
                                capture_output=True, text=True, check=True, timeout=20)
        if result.stdout.split() != [prefix, source, filesystem]:
            raise cache.Refused('Unexpected filesystem mount: ' + prefix)
        old = {r['device'] for r in rows if 'device' in r and r['path'].startswith(prefix + '/')}
        if len(old) != 1:
            raise cache.Refused('Ambiguous original filesystem identity: ' + prefix)
        mapping.append((prefix, old.pop(), os.stat(prefix).st_dev))
    return mapping


def reconcile(manifest, journal, mounts=()):
    expected = {r['path']: r for r in manifest['entries']}
    events = [json.loads(line) for line in journal.read_text().splitlines()]
    if not events or events[0].get('event') != 'start':
        raise cache.Refused('Missing original apply journal')
    intents, removed = set(), set()
    for event in events:
        kind, path = event['event'], event.get('path')
        if kind in {'intent', 'removed'}:
            if path not in expected:
                raise cache.Refused('Journal path outside approved scope')
            if kind == 'intent':
                intents.add(path)
            else:
                if path not in intents or path in removed:
                    raise cache.Refused('Invalid removal journal')
                removed.add(path)
        elif kind not in {'start', 'stopped', 'resume', 'complete'}:
            raise cache.Refused('Unrecognized journal event')
    remaining, recovered = [], []
    for path, row in expected.items():
        try:
            current = mapped_devices(dict(path=path, **cache.identity(path)), mounts)
        except FileNotFoundError:
            if path not in intents:
                raise cache.Refused('Missing entry without prior deletion intent: ' + path)
            if path not in removed:
                recovered.append(path)
            continue
        if path in removed:
            raise cache.Refused('Deleted path reappeared: ' + path)
        fields = ('device', 'inode', 'mode') if stat.S_ISDIR(row['mode']) else (
            'device', 'inode', 'mode', 'nlink', 'bytes', 'mtime_ns')
        if any(current[k] != row[k] for k in fields) or current.get('link') != row.get('link'):
            raise cache.Refused('Remaining approved entry changed: ' + path)
        remaining.append(row)
    # A no-follow enumeration detects additions, including inside nested directories.
    actual = set()
    for task in manifest['tasks']:
        if task not in {r['path'] for r in remaining}:
            continue
        if Path(task).resolve(strict=True) != Path(task):
            raise cache.Refused('Task parents became symlinks')
        stack = [task]
        while stack:
            directory = stack.pop()
            actual.add(directory)
            with os.scandir(directory) as entries:
                for entry in entries:
                    actual.add(entry.path)
                    if entry.is_dir(follow_symlinks=False):
                        stack.append(entry.path)
    if actual != {r['path'] for r in remaining}:
        raise cache.Refused('Unexpected directory contents')
    return remaining, recovered


def verify_provenance(manifest, home):
    count = 0
    for row in manifest['entries']:
        source = Path(row['path'])
        if stat.S_ISREG(row['mode']) and row['bytes'] <= cache.SMALL and source.name.startswith(('.command', '.exitcode')):
            copy = home / 'task_provenance' / source.relative_to(cache.BASE)
            if copy.is_symlink() or not copy.is_file() or copy.stat().st_size != row['bytes']:
                raise cache.Refused('Missing preserved task record: ' + str(copy))
            if row.get('sha256') and cache.sha(copy.read_bytes()) != row['sha256']:
                raise cache.Refused('Preserved task record changed: ' + str(copy))
            count += 1
    return count


def refreshed_state(manifest, mounts):
    """Refresh all consumers; reconstruct only reconciled, retired atlas rows.

    Original rows keep dependency and trace classification stable across a
    partial deletion. They never authorize extra paths: reconcile has already
    checked every remaining entry and every absent entry's journal intent.
    """
    original_scan = cache.scan
    selected = set(manifest['tasks'])
    roots = [(r, True) for r in cache.ROOTS] + [
        (cache.REPO, False), (cache.RESULTS, False),
        (cache.BASE / 'downstream', False), (cache.RESULTS.parent / 'references', False)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        pending = {(str(r), c): pool.submit(original_scan, r, c) for r, c in roots}

        def scan(root, cache_mode=False):
            key = (str(root), cache_mode)
            result = pending[key].result() if key in pending else original_scan(root, cache_mode)
            if Path(root) == cache.ROOTS[3]:
                result = dict(result)
                result['entries'] = [r for r in result['entries'] if cache.task_for(r['path']) not in selected]
                result['entries'] += mapped_devices(manifest['entries'], mounts, to_original=False)
            return result

        with patch.object(cache, 'scan', scan):
            return mapped_devices(cache.collect(), mounts)


def resume(path, digest, apply=False):
    path = Path(path)
    raw = path.read_bytes()
    if cache.sha(raw) != digest:
        raise cache.Refused('Manifest hash mismatch')
    manifest = json.loads(raw)
    if manifest['helper_sha256'] != cache.sha(Path(cache.__file__).read_bytes()):
        raise cache.Refused('Original helper changed')
    if manifest['batch'] != 'atlas_wgbs' or not cache.atlas_concurrent(manifest['review']):
        raise cache.Refused('Recovery is restricted to independent atlas retirement')
    tasks = set(manifest['tasks'])
    if not tasks or any(cache.task_for(t) != t or not t.startswith(str(cache.ROOTS[3]) + '/') for t in tasks):
        raise cache.Refused('Invalid atlas task scope')
    if any(cache.task_for(r['path']) not in tasks for r in manifest['entries']):
        raise cache.Refused('Entry outside approved task scope')
    home = path.resolve(strict=True).parent
    if any(root == home or root in home.parents for root in cache.ROOTS):
        raise cache.Refused('Journal must be outside caches')
    journal = home / 'deletions.jsonl'
    if json.loads(journal.read_text().splitlines()[0]).get('approved_sha256') != digest:
        raise cache.Refused('Journal approval does not match manifest')
    mounts = mount_mapping(manifest)
    # This recovery is executed inside the user's replacement interactive
    # allocation. Admit only that exact scheduler record, not other new jobs.
    maintenance_job = os.environ.get('SLURM_JOB_ID', '')
    maintenance_description = maintenance_job + '|vscode|/path/to/hpc-home|(null)'
    queue = cache.scheduler()
    if maintenance_job.isdigit() and maintenance_description in queue['jobs']:
        manifest['review']['concurrent_job_evidence'][cache.sha(maintenance_description.encode())] = {
            'job': maintenance_description,
            'reason': 'Current interactive allocation running approved interrupted-cleanup recovery'}
    with cache.maintenance(cache.ATLAS_LOCK):
        cache.require_idle(manifest['review'], queue)
        remaining, recovered = reconcile(manifest, journal, mounts)
        count = verify_provenance(manifest, home)
        print(json.dumps({'remaining_entries': len(remaining), 'recovered_unfinished_intents': len(recovered),
                          'preserved_records': count, 'apply': apply}), flush=True)
        if not apply:
            return
        state = refreshed_state(manifest, mounts)
        if cache.unresolved_errors(state, manifest['review']) or set(cache.batch_selection(state, 'atlas_wgbs')) != tasks:
            raise cache.Refused('Dependency scan or candidate scope changed')
        cache.check_protected(manifest['state'], state, manifest['review'])
        def findings(value):
            return sorted(cache.encode(f) for f in cache.review_findings(value, manifest['review']))
        if findings(state) != findings(manifest['state']):
            raise cache.Refused('Dependency findings changed')
        def traces(value):
            return [r for r in value['traces'] if '/bcp_all_atlas_step3/' in r['path']]
        if traces(state) != traces(manifest['state']):
            raise cache.Refused('Atlas trace evidence changed')
        errors = cache.validate_review(manifest['review'], state, manifest['tasks'])
        if errors:
            raise cache.Refused('; '.join(errors))
        remaining, recovered = reconcile(manifest, journal, mounts)
        cache.require_idle(manifest['review'], cache.scheduler())
        before = cache.allocated_usage((cache.ROOTS[3],))
        with journal.open('a') as handle:
            def log(event, **values):
                handle.write(json.dumps(dict(time=cache.now(), event=event, **values)) + '\n')
                handle.flush()
                os.fsync(handle.fileno())
            log('resume', recovery_helper_sha256=cache.sha(Path(__file__).read_bytes()),
                maintenance_allocation=maintenance_description, mount_device_mapping=mounts)
            for item in recovered:
                log('removed', path=item, reconciled_absent_after_intent=True)
            try:
                for row in sorted(remaining, key=lambda r: len(Path(r['path']).parts), reverse=True):
                    log('intent', path=row['path'])
                    cache.remove_entry(mapped_devices(row, mounts, to_original=False))
                    log('removed', path=row['path'])
                log('complete', tasks=len(tasks))
            except BaseException as exc:
                log('stopped', reason=str(exc))
                raise
        after = cache.allocated_usage((cache.ROOTS[3],))
        cache.write_new(home / 'recovery_execution.json', {
            'completed_utc': cache.now(), 'manifest_sha256': digest,
            'tasks': len(tasks), 'removed_entries_this_attempt': len(remaining),
            'reconciled_intents': len(recovered), 'allocated_before_resume': before,
            'allocated_after': after,
            'measured_reduction_this_attempt_bytes': sum(before.values()) - sum(after.values()),
            'postverification': 'pending independent retained-path and journal verification'})
        print('Approved atlas entries removed; independent verification pending.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--approved-manifest', required=True)
    parser.add_argument('--confirm-manifest-sha256', required=True)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    resume(args.approved_manifest, args.confirm_manifest_sha256, args.apply)


if __name__ == '__main__':
    main()
