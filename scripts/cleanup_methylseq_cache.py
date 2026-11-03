#!/usr/bin/env python3
"""Inventory bounded Nextflow caches; delete only exact, separately approved tasks."""
import argparse
import collections
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
BASE = Path('/path/to/hpc-scratch/tmp/em_seq')
ROOTS = tuple(BASE / x for x in ('work_emseq_integrated', 'work_emseq',
             'work_emseq_spikein', 'work_bcp_all_wgbs_step3'))
RESULTS = Path('/path/to/hpc-scratch/nextflow/nf-core-methylseq/results')
LOCK = BASE / '.methylseq-cache-maintenance.lock'
ATLAS_LOCK = BASE / '.atlas-cache-maintenance.lock'
PROTECTED_RUNS = {'2026-08-31-emseq-integrated-8samples-trim-8-3-3-10-v1',
                 '2026-09-16-emseq-integrated-heh-4new-trim-8-3-3-10-v1'}
OLD_RUNS = {'2026-07-20-emseq-integrated-v4.2.0',
            '2026-08-21-emseq-integrated-wgbstools-hg38-8samples-v1',
            '2026-08-27-emseq-integrated-8samples-zero-hard-clip-mbias-v1'}
ABANDONED = {'04/9ef79c8420f1dd3b6601ea1a9e3927', '08/9a16522dc9b9edd5c2abc25d0d757b',
             'd1/a4467768e6a24423823e76dfcf465a', 'e3/e4ba7a388e8bc39c8b8f6790c887a0'}
SMALL = 1024 * 1024
PRUNE = {'.git', '.codex', '.agents', '__pycache__', 'software', 'library', '.nextflow'}
TEXT = {'.json', '.tsv', '.csv', '.sh', '.py', '.R', '.yml', '.yaml', '.conf', '.config', '.txt', '.log'}


class Refused(ValueError):
    pass


def encode(value):
    return (json.dumps(value, sort_keys=True, indent=2) + '\n').encode()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def write_new(path, value):
    with Path(path).open('xb') as h:
        h.write(encode(value))
        h.flush()
        os.fsync(h.fileno())


def identity(path):
    s = Path(path).lstat()
    value = dict(device=s.st_dev, inode=s.st_ino, mode=s.st_mode, nlink=s.st_nlink,
                 bytes=s.st_size, allocated_bytes=s.st_blocks * 512, mtime_ns=s.st_mtime_ns)
    if stat.S_ISLNK(s.st_mode):
        value['link'] = os.readlink(path)
    return value


def inspect(root, cache=False, emit=None):
    """Read-only worker. Caller imposes a timeout on all filesystem operations."""
    root = Path(root)
    if root.is_symlink() or not root.is_dir():
        return {'root': str(root), 'entries': [], 'errors': ['Scan root is missing or a symlink: ' + str(root)]}
    rows = []
    errors = []
    stack = [root]
    while stack:
        directory = stack.pop()
        try:
            entries = sorted(os.scandir(directory), key=lambda e: e.name)
            for e in entries:
                path = Path(e.path)
                if not cache and (e.name in PRUNE or path == REPO / 'docs/audits'):
                    continue
                r = dict(path=str(path), **identity(path))
                if stat.S_ISDIR(r['mode']):
                    stack.append(path)
                elif stat.S_ISLNK(r['mode']):
                    # Resolve chains without descending into a linked directory.
                    try:
                        target = Path(os.path.normpath(os.path.join(str(path.parent), r['link'])))
                        # Raw inputs outside project-owned roots are not potential cache
                        # consumers. Do not traverse remote data-provider mounts.
                        r['link_target'] = str(target)
                        owned = (str(target).startswith('/path/to/hpc-home/') or
                                 str(target).startswith('/path/to/hpc-scratch/') or
                                 any(t == target or t in target.parents for t in ROOTS))
                        r['resolved'] = str(path.resolve(strict=True)) if owned else str(target)
                        r['target_is_dir'] = Path(r['resolved']).is_dir() if owned else False
                        if not owned:
                            r['external_target_not_traversed'] = True
                    except (OSError, RuntimeError) as exc:
                        r['link_error'] = str(exc)
                        # ENOENT outside the caches can be resolved as historical
                        # evidence by an explicit review. Permissions, loops and
                        # missing cache targets never qualify for that exception.
                        if isinstance(exc, FileNotFoundError):
                            r['missing_target'] = str(path.resolve(strict=False))
                elif stat.S_ISREG(r['mode']):
                    wanted = (cache and (e.name in {'.command.sh', '.command.run', '.exitcode'})) or (
                        not cache and path.suffix in TEXT)
                    if wanted and r['bytes'] <= SMALL:
                        data = path.read_bytes()
                        r['sha256'] = sha(data)
                        # Only dependency evidence, not arbitrary clinical text, is recorded.
                        text = data.decode(errors='replace')
                        if cache and e.name == '.command.run':
                            mounts = []
                            lines = []
                            for line in text.splitlines():
                                if 'apptainer exec ' in line or 'singularity exec ' in line:
                                    for known in ROOTS:
                                        pattern = r'(?<= )-B ' + re.escape(str(known)) + r'(?=\s)'
                                        line, count = re.subn(pattern, '-B <recorded-cache-mount>', line)
                                        if count:
                                            mounts.append(str(known))
                                lines.append(line)
                            text = '\n'.join(lines)
                            r['container_mounts'] = sorted(set(mounts))
                        r['references'] = sorted(set(re.findall(
                            re.escape(str(BASE)) + r'/[^\s\x00\"\'<>;|)]+', text)))
                        if e.name == '.exitcode':
                            r['exit'] = text.strip()
                    elif wanted:
                        r['unscanned_text'] = True
                else:
                    r['special'] = True
                rows.append(r)
                if emit:
                    emit(r)
        except OSError as exc:
            errors.append(f'{directory}: {exc}')
    return {'root': str(root), 'entries': rows, 'errors': errors}


def scan(root, cache=False):
    command = [sys.executable, __file__, '_inspect', str(root), '--stream'] + (['--cache'] if cache else [])
    errors = []
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True, timeout=900)
        output = result.stdout
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ''
        if isinstance(output, bytes):
            output = output.decode(errors='replace')
        errors.append(f'Scan timed out (partial inventory retained): {root}')
    except (OSError, subprocess.SubprocessError) as exc:
        return {'root': str(root), 'entries': [], 'errors': [f'Scan unavailable: {root}: {type(exc).__name__}']}
    entries = []
    complete = False
    for line in output.splitlines():
        try:
            row = json.loads(line)
            if '_scan_errors' in row:
                errors.extend(row['_scan_errors']); complete = True
            else:
                entries.append(row)
        except ValueError:
            errors.append(f'Incomplete scan record: {root}')
    if not complete and not errors:
        errors.append(f'Incomplete scan: {root}')
    return {'root': str(root), 'entries': entries, 'errors': errors}


def scheduler():
    try:
        p = subprocess.run(['squeue', '--noheader', '--user', str(os.getuid()),
                            '--format=%i|%j|%Z|%o'], capture_output=True, text=True,
                           check=True, timeout=20)
        return {'jobs': sorted(x for x in p.stdout.splitlines() if x.strip()), 'error': None}
    except (OSError, subprocess.SubprocessError) as exc:
        return {'jobs': [], 'error': str(exc)}


def task_for(path):
    value = str(path)
    for root in ROOTS:
        prefix = str(root) + '/'
        if not value.startswith(prefix):
            continue
        parts = value[len(prefix):].split('/', 2)
        if len(parts) >= 2 and re.fullmatch('[0-9a-f]{2}', parts[0]) and re.fullmatch('[0-9a-f]{30}', parts[1]):
            return prefix + parts[0] + '/' + parts[1]
    return None


def memberships(tasks):
    members = collections.defaultdict(set)
    errors = []
    evidence = []
    parents = (RESULTS / 'bulk', *(RESULTS / 'bcp_all_atlas_step3' / mode
                                  for mode in ('qc', 'final', 'sensitivity')))
    for parent in parents:
        for trace in sorted(parent.glob('*/pipeline_info/*trace*.txt')):
            data = trace.read_bytes()
            evidence.append({'path': str(trace), 'sha256': sha(data)})
            for row in csv.DictReader(data.decode().splitlines(), delimiter='\t'):
                short = row.get('hash', '')
                if not re.fullmatch('[0-9a-f]{2}/[0-9a-f]{6,30}', short):
                    continue
                matches = [t for t in tasks if '/'.join(Path(t).parts[-2:]).startswith(short)]
                if len(matches) > 1:
                    errors.append(f'Ambiguous trace hash {short}: {trace}')
                elif matches:
                    members[matches[0]].add(trace.parent.parent.name)
    # A killed controller can leave attempts out of its final trace. Saved
    # scheduler submission/TaskHandler records contain exact work directories.
    logfiles = sorted(logfile for parent in parents
                      for logfile in parent.glob('*/pipeline_info/nextflow_*.log'))
    for logfile in logfiles:
        digest = hashlib.sha256()
        try:
            with logfile.open('rb') as handle:
                for raw in handle:
                    digest.update(raw)
                    for match in re.finditer(rb'workDir: (\S+)', raw):
                        task = match.group(1).decode(errors='replace').rstrip(';]')
                        if task in tasks:
                            members[task].add(logfile.parent.parent.name)
            evidence.append({'path': str(logfile), 'sha256': digest.hexdigest(),
                             'kind': 'saved_task_log'})
        except OSError as exc:
            errors.append('Unreadable task membership log: ' + str(logfile) + ': ' + str(exc))
    if not PROTECTED_RUNS.issubset(set().union(*members.values()) if members else set()):
        errors.append('Missing protected-run trace membership')
    return members, errors, evidence


def ancestors_target(target, tasks):
    # Lexical prefix checks avoid repeatedly stat-ing large input trees.
    p = os.path.normpath(str(target)).rstrip('/')
    return {t for t in tasks if p == t or t.startswith(p + '/') or p.startswith(t + '/')}


def assess(cache_scans, external_scans, members):
    rows = [r for s in cache_scans for r in s['entries']]
    owned = collections.defaultdict(list)
    for r in rows:
        t = task_for(r['path'])
        if t:
            owned[t].append(r)
    tasks = set(owned)
    protected = {t for t in tasks if members.get(t, set()) & PROTECTED_RUNS}
    edges = collections.defaultdict(set)
    blockers = [e for s in cache_scans + external_scans for e in s['errors']]
    findings = []
    for r in rows + [r for s in external_scans for r in s['entries']]:
        source = task_for(r['path'])
        targets = list(r.get('references', []))
        if 'link_target' in r:
            targets.append(r['link_target'])
        if 'resolved' in r:
            targets.append(r['resolved'])
        if r.get('link_error'):
            # Never turn unresolvable links into evidence of independence.
            blockers.append('Unresolved link: ' + r['path'])
            target = r.get('missing_target')
            if not source and target and not any(
                    value == str(root) or value.startswith(str(root) + '/')
                    for value in (target, r.get('link_target', '')) for root in ROOTS):
                findings.append({'path': r['path'], 'target': target,
                                 'link': r['link'], 'kind': 'missing_external_link'})
        if r.get('external_target_not_traversed'):
            findings.append({'path': r['path'], 'target': r['resolved'], 'kind': 'external_target_not_traversed'})
        for target in r.get('container_mounts', []):
            findings.append({'path': r['path'], 'target': target,
                             'kind': 'container_mount', 'sha256': r['sha256']})
        if r.get('unscanned_text'):
            findings.append({'path': r['path'], 'kind': 'oversized_text',
                             'bytes': r['bytes'], 'mtime_ns': r['mtime_ns'],
                             'inode': r['inode'], 'device': r['device']})
        for target in targets:
            hits = ancestors_target(target, tasks)
            if source:
                edges[source].update(hits - {source})
            elif hits:
                # Provenance text needs review; resolving external links is a live dependency.
                if 'link' in r:
                    protected.update(hits)
                findings.append({'path': r['path'], 'target': target, 'tasks': sorted(hits),
                                 'kind': 'external_link' if 'link' in r else 'text_reference',
                                 **({'sha256': r['sha256']} if 'sha256' in r else {})})
    batches = {}
    for t, entries in owned.items():
        rel = '/'.join(Path(t).parts[-2:])
        root = Path(t).parents[1]
        runs = members.get(t, set())
        if root == ROOTS[0] and rel in ABANDONED and runs <= OLD_RUNS:
            exitcodes = [r.get('exit') for r in entries if Path(r['path']).name == '.exitcode']
            if exitcodes == ['143']:
                batches[t] = 'abandoned_alignment'
        elif root == ROOTS[0] and runs and runs <= OLD_RUNS:
            batches[t] = 'superseded_integrated'
        elif root in ROOTS[1:3]:
            batches[t] = 'legacy_emseq'
        elif root == ROOTS[3]:
            batches[t] = 'atlas_wgbs'
        if any(r.get('special') or (stat.S_ISREG(r['mode']) and r['nlink'] != 1) for r in entries):
            protected.add(t)
    # Unmatched and otherwise unselected tasks are retained, including their upstream inputs.
    protected.update(tasks - set(batches))
    while True:
        expanded = protected | set().union(*(edges[t] for t in protected))
        if expanded == protected:
            break
        protected = expanded
    candidates = {t: b for t, b in batches.items() if t not in protected}
    return owned, candidates, protected, blockers, findings


def collect():
    caches = []
    for root in ROOTS:
        print(f'Inspecting cache {root}', file=sys.stderr, flush=True)
        caches.append(scan(root, True))
    external = []
    # All project outputs/references, including physical targets of repository navigation links.
    pending = [REPO, RESULTS, BASE / 'downstream', RESULTS.parent / 'references']
    visited = set()
    covered = set()
    while pending:
        root = Path(pending.pop())
        if str(root) in visited or any(root == r or r in root.parents for r in ROOTS) or any(p in root.parents for p in covered):
            continue
        visited.add(str(root))
        print(f'Inspecting consumers {root}', file=sys.stderr, flush=True)
        result = scan(root)
        external.append(result)
        if not result['errors']:
            covered.add(root)
        for r in result['entries']:
            if 'resolved' in r and r['resolved'].startswith('/path/to/hpc-scratch/'):
                target = Path(r['resolved'])
                if not task_for(target) and r.get('target_is_dir') and str(target) not in visited:
                    pending.append(target)
    tasks = {task_for(r['path']) for s in caches for r in s['entries']} - {None}
    members, trace_errors, traces = memberships(tasks)
    owned, candidates, protected, errors, findings = assess(caches, external, members)
    return dict(tasks=owned, candidates=candidates, protected=sorted(protected),
                memberships={k: sorted(v) for k, v in members.items()},
                errors=errors + trace_errors, findings=findings, traces=traces,
                scan_roots=sorted(visited), queue=scheduler())


def atlas_concurrent(review):
    return review.get('maintenance_policy') == 'independent_atlas'


def review_findings(state, review):
    if not atlas_concurrent(review):
        return state['findings']
    root = str(ROOTS[3])
    # Unresolved and oversized sources still require review. Identified text or
    # mount references exclusively to other caches are irrelevant to this batch.
    return [f for f in state['findings'] if f['kind'] not in {'text_reference', 'container_mount'}
            or root in f.get('target', '') or root in f['path']
            or any(t.startswith(root + '/') for t in f.get('tasks', []))]


def check_protected(before, after, review):
    for task in before['protected']:
        old = before['tasks'][task]
        new = after['tasks'].get(task, [])
        mutable = set(before.get('memberships', {}).get(task, [])) & set(review.get('concurrent_runs', []))
        if atlas_concurrent(review) and mutable and not task.startswith(str(ROOTS[3]) + '/'):
            # Active tasks may grow outputs. Preserve their physical directories;
            # fresh full dependency scans still veto any atlas consumers.
            fields = ('device', 'inode', 'mode')
            old = [{k:r[k] for k in fields} for r in old if r['path'] == task]
            new = [{k:r[k] for k in fields} for r in new if r['path'] == task]
            if len(old) != 1 or len(new) != 1:
                raise Refused('Protected task directory missing: ' + task)
        if old != new:
            raise Refused('Protected task changed: ' + task)


def validate_review(review, state, selected):
    errors = []
    try:
        age = (datetime.now(timezone.utc) - datetime.fromisoformat(review['checked_utc'])).total_seconds()
        if not 0 <= age <= 86400:
            errors.append('Dependency review must be within the last 24 hours')
        if set(review['tasks']) != set(selected) or review['recovery_policy'] != 'retire_resume':
            errors.append('Review task scope or retirement policy mismatch')
        if atlas_concurrent(review) and any(not t.startswith(str(ROOTS[3]) + '/') for t in selected):
            errors.append('Independent maintenance is restricted to the atlas cache')
        if atlas_concurrent(review) and not set(review.get('concurrent_runs', [])) <= PROTECTED_RUNS:
            errors.append('Unrecognized concurrent protected run')
        for task in selected:
            item = review['tasks'][task]
            if item.get('no_active_consumers') is not True or not item.get('evidence'):
                errors.append('Missing task-specific release evidence: ' + task)
        for f in review_findings(state, review):
            fid = sha(encode(f))
            entry = review['findings'].get(fid, {})
            if entry.get('status') not in {'historical', 'released'} or not entry.get('reason'):
                errors.append('Unreviewed dependency finding: ' + fid)
        evidence = review['retained_evidence']
        if not evidence:
            errors.append('Missing retained-output/reference evidence')
        for item in evidence:
            p = Path(item['path'])
            if not p.is_file() or p.is_symlink() or p.stat().st_size > 32 * SMALL or sha(p.read_bytes()) != item['sha256']:
                errors.append('Retained evidence changed: ' + str(p))
        # Atlas output completeness must be explicitly supported, not inferred from absence of links.
        if any(state['candidates'][t] == 'atlas_wgbs' for t in selected):
            if review.get('atlas_published_inputs_independent') is not True:
                errors.append('Atlas published-input independence not established')
    except (KeyError, TypeError, ValueError, OSError) as exc:
        errors.append('Invalid review: ' + str(exc))
    return errors


def unresolved_errors(state, review):
    """Only an individually reviewed, missing non-cache consumer link is exempt."""
    released = set()
    for finding in state['findings']:
        if finding['kind'] != 'missing_external_link':
            continue
        decision = (review or {}).get('findings', {}).get(sha(encode(finding)), {})
        if decision.get('status') == 'historical' and decision.get('reason'):
            released.add('Unresolved link: ' + finding['path'])
    return [e for e in state['errors'] if e not in released]


def metadata_destination(path):
    path = Path(path).absolute()
    if '..' in path.parts or path.resolve(strict=False) != path:
        raise Refused('Metadata destination must have physical, unaliased parents')
    if path.exists() or path.is_symlink():
        raise Refused('Metadata destination must be new')
    scratch = Path('/path/to/hpc-scratch/tmp')
    if scratch not in path.parents or any(root == path or root in path.parents for root in ROOTS):
        raise Refused('Put audit metadata in scratch tmp outside the work caches')
    return path


def batch_selection(state, batch):
    selected = {t for t, b in state['candidates'].items() if b == batch}
    # Other candidate batches remain retained until their own approval/apply.
    while True:
        hits = set()
        for task, entries in state['tasks'].items():
            if task in selected:
                continue
            for r in entries:
                targets = r.get('references', []) + [r[k] for k in ('link_target', 'resolved') if k in r]
                for target in targets:
                    hits.update(ancestors_target(target, selected))
        reduced = selected - hits
        if reduced == selected:
            return sorted(selected)
        selected = reduced


def plan(destination, batch, review_path=None, state=None):
    destination = metadata_destination(destination)
    state = collect() if state is None else state
    selected = batch_selection(state, batch)
    review = json.loads(Path(review_path).read_text()) if review_path else None
    blockers = unresolved_errors(state, review)
    if not selected:
        blockers.append('No eligible task directories')
    if state['queue']['error']:
        blockers.append('Scheduler unavailable')
    blockers += validate_review(review, state, selected) if review else ['No reviewed dependency release']
    records = [r for t in selected for r in state['tasks'][t]]
    manifest = dict(schema=1, created_utc=now(), batch=batch, tasks=selected, entries=records,
                    state=state, review=review, blockers=sorted(set(blockers)),
                    helper_sha256=sha(Path(__file__).read_bytes()))
    destination.mkdir(mode=0o700, parents=True)
    write_new(destination / 'manifest.json', manifest)
    digest = sha((destination / 'manifest.json').read_bytes())
    (destination / 'manifest.sha256').write_text(digest + '  manifest.json\n')
    write_new(destination / 'summary.json', dict(tasks=len(selected), entries=len(records),
              allocated_bytes=sum(r['allocated_bytes'] for r in records if stat.S_ISREG(r['mode'])),
              manifest_sha256=digest, blockers=manifest['blockers'], apply_requires_idle=True))
    print(json.dumps({'manifest': str(destination / 'manifest.json'), 'sha256': digest,
                      'tasks': len(selected), 'blockers': len(manifest['blockers'])}))


@contextmanager
def maintenance(lock=None):
    with (LOCK if lock is None else lock).open('a') as h:
        fcntl.flock(h, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def require_idle(review, queue):
    if queue['error']:
        raise Refused('Cannot establish scheduler clearance')
    allowed = review.get('unrelated_jobs', {})
    for job in queue['jobs']:
        if atlas_concurrent(review):
            if any(s in job.lower() for s in ('work_bcp_all_wgbs_step3', 'bcp_all', 'bcpall', 'atlas', 'wgbs')):
                raise Refused('Atlas-related queued/running job blocks apply: ' + job)
            if not review.get('concurrent_job_evidence', {}).get(sha(job.encode())):
                raise Refused('Unreviewed concurrent job: ' + job)
            continue
        # Explicit descriptions bind unrelated allocations, never pipeline/work consumers.
        if any(s in job.lower() for s in ('nextflow', 'methylseq', 'bismark', 'dmrichr', 'bcp_all', 'emseq_pat', 'emseq_uxm', 'emseq_epidish')):
            raise Refused('Relevant queued/running job blocks cache apply: ' + job)
        if not allowed.get(sha(job.encode())):
            raise Refused('Unreviewed current job: ' + job)


def verify_entries(entries, tasks):
    actual = []
    for task in tasks:
        if task_for(task) != task or Path(task).resolve(strict=True) != Path(task):
            raise Refused('Invalid physical task path')
        result = scan(task, True)
        if result['errors']:
            raise Refused('Cannot re-inventory task')
        actual += [dict(path=task, **identity(task))] + result['entries']
    # Application-independent identity fields; dependency metadata is rechecked separately.
    def basic(row):
        return {k: row[k] for k in ('path', 'device', 'inode', 'mode', 'nlink', 'bytes', 'mtime_ns', 'link') if k in row}
    if sorted(map(basic, actual), key=lambda r:r['path']) != sorted(map(basic, entries), key=lambda r:r['path']):
        raise Refused('Task contents or identities changed')


@contextmanager
def open_parent(path):
    fd = os.open('/', os.O_RDONLY | os.O_DIRECTORY)
    try:
        for part in Path(path).parent.parts[1:]:
            child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
            os.close(fd)
            fd = child
        yield fd
    finally:
        os.close(fd)


def remove_entry(row):
    p = Path(row['path'])
    with open_parent(p) as fd:
        st = os.stat(p.name, dir_fd=fd, follow_symlinks=False)
        if (st.st_dev, st.st_ino, st.st_mode) != (row['device'], row['inode'], row['mode']):
            raise Refused('Entry replaced before removal: ' + str(p))
        if stat.S_ISDIR(st.st_mode):
            os.rmdir(p.name, dir_fd=fd)  # Refuses unexpected children.
        else:
            if (st.st_size, st.st_mtime_ns, st.st_nlink) != (row['bytes'], row['mtime_ns'], row['nlink']):
                raise Refused('Entry changed before removal: ' + str(p))
            if stat.S_ISLNK(st.st_mode):
                if os.readlink(p.name, dir_fd=fd) != row['link']:
                    raise Refused('Symlink changed')
            elif not stat.S_ISREG(st.st_mode) or st.st_nlink != 1:
                raise Refused('Nonregular or shared file refused')
            os.unlink(p.name, dir_fd=fd)
        os.fsync(fd)


def allocated_usage(roots=None):
    result = subprocess.run(['du', '-s', '-B1', '--', *map(str, ROOTS if roots is None else roots)],
                            check=True, capture_output=True, text=True, timeout=180)
    return {p: int(b) for b, p in (line.split('\t', 1) for line in result.stdout.splitlines())}


def apply(manifest_path, approved_hash):
    raw = Path(manifest_path).read_bytes()
    if sha(raw) != approved_hash:
        raise Refused('Approved manifest hash mismatch')
    m = json.loads(raw)
    if m['helper_sha256'] != sha(Path(__file__).read_bytes()) or m['blockers']:
        raise Refused('Blocked manifest or changed helper')
    if not m['tasks'] or any(task_for(t) != t for t in m['tasks']):
        raise Refused('Invalid task scope')
    home = Path(manifest_path).resolve(strict=True).parent
    if any(home == root or root in home.parents for root in ROOTS):
        raise Refused('Approval and journal must be outside all cache roots')
    concurrent = atlas_concurrent(m['review'])
    if concurrent and (m['batch'] != 'atlas_wgbs' or any(
            not t.startswith(str(ROOTS[3]) + '/') for t in m['tasks'])):
        raise Refused('Independent maintenance is restricted to the atlas batch')
    with maintenance(ATLAS_LOCK if concurrent else None):
        require_idle(m['review'], scheduler())
        state = collect()
        if unresolved_errors(state, m['review']) or set(m['tasks']) != set(batch_selection(state, m['batch'])):
            raise Refused('Candidates or scan coverage changed')
        # Stable protected identities are checked in full, even though no protected files are deleted.
        check_protected(m['state'], state, m['review'])
        def relevant_traces(s):
            return [r for r in s['traces'] if not concurrent or '/bcp_all_atlas_step3/' in r['path']]
        if ((not concurrent and state['protected'] != m['state']['protected']) or
                review_findings(state, m['review']) != review_findings(m['state'], m['review']) or
                relevant_traces(state) != relevant_traces(m['state'])):
            raise Refused('Dependency evidence changed')
        errors = validate_review(m['review'], state, m['tasks'])
        if errors:
            raise Refused('; '.join(errors))
        verify_entries(m['entries'], m['tasks'])
        require_idle(m['review'], scheduler())
        measure_roots = (ROOTS[3],) if concurrent else None
        before_usage = allocated_usage(measure_roots)
        provenance = home / 'task_provenance'
        provenance.mkdir(mode=0o700, exist_ok=False)
        for r in m['entries']:
            p = Path(r['path'])
            if stat.S_ISREG(r['mode']) and r['bytes'] <= SMALL and p.name.startswith(('.command', '.exitcode')):
                target = provenance / p.relative_to(BASE)
                target.parent.mkdir(parents=True, exist_ok=True)
                data = p.read_bytes()
                target.write_bytes(data)
                current = identity(p)
                if sha(target.read_bytes()) != sha(data) or current != {k:r[k] for k in current}:
                    raise Refused('Provenance copy verification failed')
        journal = home / 'deletions.jsonl'
        with journal.open('x') as h:
            def log(event, **values):
                h.write(json.dumps(dict(time=now(), event=event, **values)) + '\n')
                h.flush(); os.fsync(h.fileno())
            log('start', approved_sha256=approved_hash)
            try:
                for row in sorted(m['entries'], key=lambda r:len(Path(r['path']).parts), reverse=True):
                    log('intent', path=row['path'])
                    remove_entry(row)
                    log('removed', path=row['path'])
                if any(Path(t).exists() or Path(t).is_symlink() for t in m['tasks']):
                    raise Refused('Retired task remains present')
                log('complete', tasks=len(m['tasks']))
            except BaseException as exc:
                log('stopped', reason=str(exc))
                raise
        post = collect()
        failures = unresolved_errors(post, m['review'])
        try:
            check_protected(state, post, m['review'])
        except Refused as exc:
            failures.append(str(exc))
        def live_links(findings):
            return {(f['path'], f.get('target')) for f in findings
                    if f['kind'] == 'external_link' and task_for(f['path']) not in m['tasks']
                    and (not concurrent or str(ROOTS[3]) in f.get('target', ''))}
        if live_links(post['findings']) != live_links(state['findings']):
            failures.append('Live consumer links changed after apply')
        if failures:
            write_new(home / 'postcheck_failed.json', failures)
            raise Refused('Post-apply checks failed; inspect postcheck_failed.json')
        after_usage = allocated_usage(measure_roots)
        write_new(home / 'execution.json', dict(completed_utc=now(), tasks=m['tasks'],
                  removed_entries=len(m['entries']), allocated_before=before_usage, allocated_after=after_usage,
                  measured_reduction_bytes=sum(before_usage.values())-sum(after_usage.values()),
                  protected_verification=('protected directory identities and atlas dependencies matched; unrelated active outputs may change'
                                          if concurrent else 'pre/post-apply matched')))


def main():
    if len(sys.argv)>1 and sys.argv[1] == '_inspect':
        if '--stream' in sys.argv[3:]:
            result = inspect(sys.argv[2], '--cache' in sys.argv[3:],
                             emit=lambda row: print(json.dumps(row), flush=True))
            print(json.dumps({'_scan_errors': result['errors']}), flush=True)
        else:
            print(json.dumps(inspect(sys.argv[2], '--cache' in sys.argv[3:])))
        return
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest-dir')
    p.add_argument('--batch', choices=['all', 'abandoned_alignment', 'superseded_integrated', 'legacy_emseq', 'atlas_wgbs'])
    p.add_argument('--dependency-review')
    p.add_argument('--apply', action='store_true')
    p.add_argument('--approved-manifest')
    p.add_argument('--confirm-manifest-sha256')
    a = p.parse_args()
    if a.apply:
        if not a.approved_manifest or not a.confirm_manifest_sha256 or a.manifest_dir or a.batch or a.dependency_review:
            p.error('Apply requires only an approved manifest and its exact hash')
        apply(a.approved_manifest, a.confirm_manifest_sha256)
    else:
        if not a.manifest_dir or not a.batch or a.approved_manifest or a.confirm_manifest_sha256:
            p.error('Dry-run requires manifest-dir and batch; approval flags require apply')
        if a.batch == 'all':
            if a.dependency_review:
                p.error('Reviewed manifests require one explicit batch')
            destination = metadata_destination(a.manifest_dir)
            state = collect()
            destination.mkdir(mode=0o700, parents=True)
            for batch in ('abandoned_alignment', 'superseded_integrated', 'legacy_emseq', 'atlas_wgbs'):
                plan(destination / batch, batch, state=state)
        else:
            plan(a.manifest_dir, a.batch, a.dependency_review)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        sys.exit('Cache cleanup refused: ' + str(exc))
