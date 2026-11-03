#!/usr/bin/env python3
"""Plan bounded published-output cleanup; apply only an explicitly approved manifest.

No cache, raw-data, reference, PAT or whole-directory deletion is supported.
See docs/audits/methylseq_cleanup_usage.md for dependency-review requirements.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import io
import json
import os
import pwd
import re
import stat
import subprocess
from contextlib import ExitStack, contextmanager
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RESULTS = Path('/path/to/hpc-scratch/nextflow/nf-core-methylseq/results')
BULK = RESULTS / 'bulk'
SCRATCH_METADATA = Path('/path/to/hpc-scratch/tmp')
JULY = {'2026-07-03-emseq-v4.2.0', '2026-07-03-emseq-v4.2.0_spikein'}
E27 = '2026-08-27-emseq-integrated-8samples-zero-hard-clip-mbias-v1'
E20_E21 = {'2026-07-20-emseq-integrated-v4.2.0',
           '2026-08-21-emseq-integrated-wgbstools-hg38-8samples-v1'}
ALLOWED_RUNS = JULY | {E27} | E20_E21
# These constants have no CLI or environment override. Tests patch them in-process.
SCAN_ROOTS = (REPO, Path('/path/to/hpc-scratch/tmp/em_seq/downstream'),
              RESULTS / 'uxm_deconvolution', RESULTS / 'epidishr', RESULTS / 'pcdh_epiclone')
SMALL_FILE_LIMIT = 16 * 1024 * 1024
SCAN_LIMIT = 2 * 1024 * 1024
TEXT_SUFFIXES = {'.sh', '.py', '.R', '.conf', '.json', '.tsv', '.csv', '.yaml', '.yml', '.md', '.log', '.txt'}
PRUNE = {'.git', '.codex', '.agents', '.nextflow', 'software', 'references',
         'RData', 'library', 'R', 'fastq', 'fastqs', 'pat', '__pycache__'}
FIELDS = ['run', 'path', 'artifact', 'device', 'inode', 'nlink', 'mode', 'bytes',
          'allocated_bytes', 'mtime_ns', 'parent_device', 'parent_inode']


class CleanupError(RuntimeError):
    """A required safety invariant is not established."""


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def digest_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def json_bytes(value) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + '\n').encode()


def safe_text(value: str) -> str:
    if any(c in value for c in ('\n', '\r', '\t', '\0', '\\')):
        raise CleanupError(f'Unsupported control/escape character in path: {value!r}')
    return value


def identity(path: Path) -> dict:
    info = path.lstat()
    parent = path.parent.stat()
    return dict(device=info.st_dev, inode=info.st_ino, nlink=info.st_nlink,
                mode=info.st_mode, bytes=info.st_size, allocated_bytes=info.st_blocks * 512,
                mtime_ns=info.st_mtime_ns, parent_device=parent.st_dev, parent_inode=parent.st_ino)


def hash_file(path: Path) -> str:
    before = identity(path)
    if not stat.S_ISREG(before['mode']):
        raise CleanupError(f'Checksum input is not a regular file: {path}')
    value = hashlib.sha256()
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    with os.fdopen(fd, 'rb') as handle:
        opened = os.fstat(handle.fileno())
        if (opened.st_dev, opened.st_ino) != (before['device'], before['inode']):
            raise CleanupError(f'Checksum input changed: {path}')
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(chunk)
    if identity(path) != before:
        raise CleanupError(f'File changed while hashing: {path}')
    return value.hexdigest()


def physical(path: str) -> Path:
    safe_text(path)
    if not path or not os.path.isabs(path) or '..' in path.split('/') or re.search(r'[*?\[\]]', path):
        raise CleanupError(f'Require an explicit absolute physical path: {path!r}')
    candidate = Path(path)
    if not candidate.exists() or candidate.resolve(strict=True) != candidate or candidate.is_symlink():
        raise CleanupError(f'Unresolved or symlink path refused: {path}')
    return candidate


def validate_runs(arguments: list[str]) -> list[Path]:
    runs = sorted(set(physical(arg) for arg in arguments))
    if not runs:
        raise CleanupError('At least one explicit --run-dir is required')
    for run in runs:
        if not run.is_dir() or run.parent != BULK or run.name not in ALLOWED_RUNS:
            raise CleanupError(f'Run outside finite allowed scope (E31 and caches are protected): {run}')
    return runs


def classes(value: str) -> list[str]:
    selected = value.split(',')
    if not selected or set(selected) - {'bam', 'context-calls'} or len(set(selected)) != len(selected):
        raise CleanupError('--artifacts must be bam, context-calls, or bam,context-calls')
    return sorted(selected)


def walk_error(error):
    raise CleanupError(f'Cannot inspect directory: {error}')


def walk_files(root: Path):
    """Never descend into symlink directories; return them as retained entries."""
    for directory, dirs, files in os.walk(root, followlinks=False, onerror=walk_error):
        dirs.sort(); files.sort()
        for name in list(dirs):
            path = Path(directory) / name
            if path.is_symlink():
                dirs.remove(name)
                yield path
        for name in files:
            yield Path(directory) / name


def artifact_class(relative: str, records: dict[str, dict], selected: list[str]) -> str | None:
    path = Path(relative)
    record = records[relative]
    if any(word in path.name.lower() for word in ('report', 'mbias', 'coverage', 'bedgraph', 'reference', 'genome', 'fastq', 'lbeta', 'cache')):
        return None
    if not stat.S_ISREG(record['mode']) or record['bytes'] == 0:
        return None
    if ('bam' in selected and str(path.parent) in {'bismark/alignments', 'bismark/deduplicated'}
            and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*\.bam(?:\.bai)?', path.name)):
        partner = relative[:-4] if relative.endswith('.bam.bai') else relative + '.bai'
        other = records.get(partner)
        if other and stat.S_ISREG(other['mode']) and other['bytes'] > 0:
            return 'bam'
    if ('context-calls' in selected and str(path.parent) == 'bismark/methylation_calls/methylation_calls'
            and re.fullmatch(r'(?:CpG|CHG|CHH)_[A-Za-z0-9][A-Za-z0-9._-]*\.txt\.gz', path.name)):
        return 'context-calls'
    return None


def inventory(runs: list[Path], selected: list[str]) -> tuple[list[dict], list[dict], list[str]]:
    candidates, retained, blockers = [], [], []
    for run in runs:
        records = {safe_text(str(p.relative_to(run))): identity(p) for p in walk_files(run)}
        if not any(name.startswith('pipeline_info/') and stat.S_ISREG(row['mode']) and row['bytes'] > 0
                   for name, row in records.items()):
            blockers.append(f'{run.name}: missing nonempty pipeline_info provenance')
        for relative, record in sorted(records.items()):
            kind = artifact_class(relative, records, selected)
            row = {'run': str(run), 'path': relative, 'artifact': kind or 'retained', **record}
            path = run / relative
            if kind:
                if record['nlink'] != 1:
                    raise CleanupError(f'Selected file has multiple hard links: {path}')
                if path.resolve(strict=True) != path:
                    raise CleanupError(f'Selected path has a symlink component: {path}')
                candidates.append(row)
            else:
                if stat.S_ISLNK(record['mode']):
                    row['link_target'] = safe_text(os.readlink(path))
                elif stat.S_ISREG(record['mode']) and record['bytes'] <= SMALL_FILE_LIMIT:
                    row['sha256'] = hash_file(path)
                retained.append(row)
    return candidates, retained, blockers


def scan_dependencies(runs: list[Path]) -> list[dict]:
    """A bounded evidence search, never a claim that external consumers are absent."""
    findings = {}
    for root in SCAN_ROOTS:
        if not root.exists():
            continue
        for directory, dirs, files in os.walk(root, followlinks=False, onerror=walk_error):
            directory = Path(directory)
            for name in sorted(list(dirs)):
                path = directory / name
                if path.is_symlink():
                    files.append(name); dirs.remove(name)
                elif name in PRUNE or name.startswith('runtime-') or path == REPO / 'docs/audits':
                    dirs.remove(name)
            for name in sorted(files):
                path = directory / name
                kind, matches = '', []
                if path.is_symlink():
                    resolved = path.resolve()
                    matches = [r.name for r in runs if resolved == r or r in resolved.parents]
                    kind = 'symlink'
                    fingerprint = digest_bytes(os.readlink(path).encode())
                elif path.suffix in TEXT_SUFFIXES and path.is_file():
                    if path.stat().st_size > SCAN_LIMIT:
                        kind, matches = 'unscanned-large-text', [r.name for r in runs]
                        fingerprint = digest_bytes(json_bytes(identity(path)))
                    else:
                        data = path.read_bytes()
                        matches = [r.name for r in runs if r.name.encode() in data or str(r).encode() in data]
                        fingerprint = digest_bytes(data)
                        kind = 'text-reference'
                if matches:
                    item = {'path': str(path), 'kind': kind, 'runs': sorted(matches), 'sha256': fingerprint}
                    item['id'] = digest_bytes(json_bytes(item))
                    findings[item['id']] = item
    return sorted(findings.values(), key=lambda row: row['id'])


def scheduler() -> dict:
    # No elapsed-time field: unchanged unrelated jobs keep stable fingerprints.
    command = ['squeue', '--noheader', '--user', pwd.getpwuid(os.getuid()).pw_name,
               '--format=%i|%j|%Z|%o']
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=20)
        lines = sorted(line.strip() for line in result.stdout.splitlines() if line.strip())
        return {'checked_utc': now(), 'jobs': [{'id': line.split('|', 1)[0],
                  'sha256': digest_bytes(line.encode()), 'description': line} for line in lines], 'error': None}
    except (OSError, subprocess.SubprocessError) as exc:
        return {'checked_utc': now(), 'jobs': [], 'error': str(exc)}


def read_small_json(path: Path) -> dict:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 32 * 1024 * 1024:
        raise CleanupError(f'Invalid or oversized metadata file: {path}')
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise CleanupError(f'Expected a JSON object: {path}')
    return value


def review_blockers(review: dict | None, runs: list[Path], findings: list[dict], queue: dict) -> list[str]:
    reasons = []
    if not review:
        return ['No reviewed dependency release supplied; inventory is not deletion approval']
    try:
        checked = datetime.fromisoformat(review['checked_utc'])
        age = (datetime.now(timezone.utc) - checked).total_seconds()
        if not 0 <= age <= 86400:
            reasons.append('Dependency review must be dated within the last 24 hours')
        release = review['runs']
        if set(release) != {str(r) for r in runs}:
            reasons.append('Dependency review run set differs from selected runs')
        for run in runs:
            decision = release.get(str(run), {})
            if decision.get('no_active_consumers') is not True or not decision.get('evidence'):
                reasons.append(f'{run.name}: no documented active-consumer release')
            if decision.get('recovery_policy') not in {'accept_no_guaranteed_replay', 'verified_backup'}:
                reasons.append(f'{run.name}: recovery/retention policy not reviewed')
            if decision.get('recovery_policy') == 'verified_backup' and not decision.get('backup_evidence'):
                reasons.append(f'{run.name}: verified backup evidence missing')
            if run.name == E27 and decision.get('mbias_closed') is not True:
                reasons.append('E27 M-bias/marker comparison closure is required')
            if run.name in E20_E21:
                for gate in ('e31_downstream_migration_verified', 'dmrichr_local_closure_verified',
                             'retained_consumers_independent_of_candidates'):
                    if decision.get(gate) is not True:
                        reasons.append(f'{run.name}: {gate} release is required')
            metadata = decision.get('historical_metadata', {})
            required = {'pipeline_version', 'params', 'config', 'samplesheet', 'reference',
                        'scientific_settings', 'resume', 'missing_originals'}
            if not required.issubset(metadata):
                reasons.append(f'{run.name}: historical provenance/limitations not recorded')
            external = decision.get('external_records', [])
            if not external:
                reasons.append(f'{run.name}: preserved external provenance evidence missing')
            if (decision.get('recovery_policy') == 'verified_backup'
                    and decision.get('backup_evidence') not in [r.get('path') for r in external]):
                reasons.append(f'{run.name}: backup evidence must be a checksummed external record')
            for record in external:
                path = physical(record['path'])
                if path.stat().st_size > SMALL_FILE_LIMIT or hash_file(path) != record['sha256']:
                    reasons.append(f'External retained evidence changed or too large: {path}')
        classified = review.get('findings', {})
        for finding in findings:
            decision = classified.get(finding['id'], {})
            if decision.get('status') not in {'historical', 'released'} or not decision.get('reason'):
                reasons.append(f'Unresolved consumer evidence: {finding["id"]} {finding["path"]}')
        if queue['error']:
            reasons.append(f'Cannot establish scheduler state: {queue["error"]}')
        for job in queue['jobs']:
            decision = review.get('unrelated_jobs', {}).get(job['id'], {})
            if decision.get('sha256') != job['sha256'] or not decision.get('reason'):
                reasons.append(f'Unreviewed current Slurm job: {job["id"]}')
    except (KeyError, TypeError, ValueError, OSError, CleanupError) as exc:
        reasons.append(f'Invalid dependency review: {exc}')
    return reasons


def tsv_bytes(rows: list[dict], fields: list[str]) -> bytes:
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    return output.getvalue().encode()


def write_new(path: Path, data: bytes):
    with path.open('xb') as handle:
        handle.write(data); handle.flush(); os.fsync(handle.fileno())
    with open_directory(path.parent) as directory_fd:
        os.fsync(directory_fd)


def totals(rows: list[dict]) -> dict:
    unique = {(row['device'], row['inode']): row for row in rows}
    return {'files': len(rows), 'unique_inodes': len(unique),
            'logical_bytes': sum(row['bytes'] for row in unique.values()),
            'allocated_bytes': sum(row['allocated_bytes'] for row in unique.values())}


def validate_metadata_destination(destination: Path):
    if destination == RESULTS or RESULTS in destination.parents:
        raise CleanupError('Manifest output must be outside scientific results')
    if REPO / 'docs/audits' in destination.parents:
        return
    if SCRATCH_METADATA in destination.parents:
        first = destination.relative_to(SCRATCH_METADATA).parts[0]
        if first.startswith(('cleanup_', 'emseq_cleanup_', 'emseq-cleanup-test-')):
            return
    raise CleanupError('Manifest metadata must be under repository docs/audits or a separate scratch tmp/cleanup_* directory')


def make_plan(args):
    runs, selected = validate_runs(args.run_dir), classes(args.artifacts)
    destination = Path(args.manifest_dir).resolve()
    validate_metadata_destination(destination)
    if destination.exists():
        raise CleanupError(f'Refusing to overwrite an earlier plan: {destination}')
    candidates, retained, inventory_blockers = inventory(runs, selected)
    blockers = list(inventory_blockers)
    findings, queue = scan_dependencies(runs), scheduler()
    review = read_small_json(Path(args.dependency_review)) if args.dependency_review else None
    blockers += review_blockers(review, runs, findings, queue)
    if not candidates:
        blockers.append('No eligible artifact files selected')
    retained_sums = ''.join(f'{r["sha256"]}  {r["run"]}/{r["path"]}\n' for r in retained if 'sha256' in r).encode()
    dependency_rows = [{'id': f['id'], 'kind': f['kind'], 'path': safe_text(f['path']),
                        'runs': ','.join(f['runs']), 'sha256': f['sha256']} for f in findings]
    deps = tsv_bytes(dependency_rows, ['id', 'kind', 'path', 'runs', 'sha256'])
    context = {'schema': 1, 'created_utc': now(), 'runs': [str(r) for r in runs],
               'artifacts': selected, 'helper_sha256': hash_file(Path(__file__)),
               'retained': retained, 'findings': findings, 'scheduler': queue,
               'dependency_review': review, 'blockers': blockers,
               'retained_sha256': digest_bytes(retained_sums), 'dependencies_sha256': digest_bytes(deps),
               'before_totals': totals(candidates + retained),
               'scan_scope': [str(root) for root in SCAN_ROOTS],
               'limitations': 'Bounded text/symlink scan; human release must cover unscanned binary, archived, remote and scheduler-external consumers. Missing historical originals are not reconstructed.'}
    context_data = json_bytes(context)
    plan = ('# context_sha256\t' + digest_bytes(context_data) + '\n').encode() + tsv_bytes(candidates, FIELDS)
    summary = {'mode': 'dry-run', 'plan_sha256': digest_bytes(plan), 'selected': totals(candidates),
               'retained': totals(retained), 'blocked': blockers,
               'note': 'Candidate bytes, not measured quota reclamation. No data changed. Apply requires explicit user approval of this exact plan hash.'}
    # Detect source changes during the scan before publishing any plan.
    if inventory(runs, selected) != (candidates, retained, inventory_blockers):
        raise CleanupError('Run contents changed while planning')
    destination.mkdir(parents=True)
    for name, data in [('context.json', context_data), ('plan.tsv', plan), ('retained.sha256', retained_sums),
                       ('dependencies.tsv', deps), ('summary.json', json_bytes(summary))]:
        write_new(destination / name, data)
    print(json.dumps(summary, indent=2))


@contextmanager
def open_directory(path: Path):
    """Open every path component without following symlinks."""
    fd = os.open('/', os.O_RDONLY | os.O_DIRECTORY)
    try:
        for part in path.parts[1:]:
            next_fd = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
            os.close(fd); fd = next_fd
        yield fd
    finally:
        os.close(fd)


def unlink_checked(row: dict):
    path = Path(row['run']) / row['path']
    with open_directory(path.parent) as parent_fd:
        parent = os.fstat(parent_fd)
        if (parent.st_dev, parent.st_ino) != (row['parent_device'], row['parent_inode']):
            raise CleanupError(f'Parent directory changed: {path}')
        info = os.stat(path.name, dir_fd=parent_fd, follow_symlinks=False)
        actual = dict(device=info.st_dev, inode=info.st_ino, nlink=info.st_nlink,
                      mode=info.st_mode, bytes=info.st_size, allocated_bytes=info.st_blocks * 512,
                      mtime_ns=info.st_mtime_ns)
        if any(actual[key] != row[key] for key in actual) or not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise CleanupError(f'Selected file changed immediately before unlink: {path}')
        os.unlink(path.name, dir_fd=parent_fd)


def sync_unlink_parent(path: Path):
    with open_directory(path.parent) as parent_fd:
        os.fsync(parent_fd)


def apply_plan(args):
    runs, selected = validate_runs(args.run_dir), classes(args.artifacts)
    plan_path = physical(args.approved_manifest)
    if plan_path.name != 'plan.tsv' or not re.fullmatch(r'[0-9a-f]{64}', args.confirm_manifest_sha256):
        raise CleanupError('Apply requires plan.tsv and its exact lowercase SHA-256')
    if plan_path.stat().st_size > 32 * 1024 * 1024:
        raise CleanupError('Manifest exceeds the metadata size limit')
    plan_data = plan_path.read_bytes()
    if digest_bytes(plan_data) != args.confirm_manifest_sha256:
        raise CleanupError('Approved manifest hash mismatch')
    directory = plan_path.parent
    validate_metadata_destination(directory)
    context_path = directory / 'context.json'
    context = read_small_json(context_path)
    header, body = plan_data.split(b'\n', 1)
    if header != ('# context_sha256\t' + hash_file(context_path)).encode():
        raise CleanupError('Bound retention/dependency context changed')
    if (context.get('schema') != 1 or context['runs'] != [str(r) for r in runs]
            or context['artifacts'] != selected or context['helper_sha256'] != hash_file(Path(__file__))
            or context['scan_scope'] != [str(root) for root in SCAN_ROOTS]):
        raise CleanupError('Helper, run scope, artifact classes or dependency scan scope changed')
    for name, key in [('retained.sha256', 'retained_sha256'), ('dependencies.tsv', 'dependencies_sha256')]:
        if hash_file(directory / name) != context[key]:
            raise CleanupError(f'Bound companion record changed: {name}')
    if context['blockers']:
        raise CleanupError('Manifest contains unresolved blockers; create and approve a new reviewed plan')
    for name in ('deletions.jsonl', 'execution_summary.json'):
        if os.path.lexists(directory / name):
            raise CleanupError('An apply attempt already exists; manual reconciliation is required')
    for run in runs:
        if os.path.lexists(run / 'PAYLOAD_RETENTION.md'):
            raise CleanupError(f'Existing retention note needs separate review: {run}')
    with ExitStack() as stack:
        # Directory locks serialize this helper only, not scientific writers.
        for run in runs:
            fd = stack.enter_context(open_directory(run))
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        candidates, retained, blockers = inventory(runs, selected)
        if tsv_bytes(candidates, FIELDS) != body or retained != context['retained'] or blockers:
            raise CleanupError('Selected or retained run contents changed since approval')
        findings = scan_dependencies(runs)
        if findings != context['findings']:
            raise CleanupError('Dependency evidence changed since approval')
        reasons = review_blockers(context['dependency_review'], runs, findings, scheduler())
        if reasons:
            raise CleanupError('\n'.join(reasons))
        removed, failure, retained_verified = [], None, False
        with (directory / 'deletions.jsonl').open('x') as journal:
            with open_directory(directory) as journal_directory_fd:
                os.fsync(journal_directory_fd)
            def log(event, row, detail=None):
                journal.write(json.dumps({'utc': now(), 'event': event, 'run': row['run'],
                                          'path': row['path'], 'detail': detail}) + '\n')
                journal.flush(); os.fsync(journal.fileno())
            try:
                for row in candidates:
                    log('intent', row)
                    try:
                        unlink_checked(row)
                    except BaseException as exc:
                        log('failed', row, str(exc)); raise
                    removed.append(row)
                    log('unlinked', row)
                    try:
                        sync_unlink_parent(Path(row['run']) / row['path'])
                    except BaseException as exc:
                        log('durability-error', row, str(exc)); raise
                _, after_retained, _ = inventory(runs, selected)
                if after_retained != retained:
                    raise CleanupError('Retained-file verification failed after deletion')
                retained_verified = True
            except BaseException as exc:
                failure = f'{type(exc).__name__}: {exc}'
            finally:
                # This manifest is single-use even if interrupted or partially applied.
                after_rows = []
                for run in runs:
                    after_rows.extend(identity(p) for p in walk_files(run))
                result = {'utc': now(), 'status': 'failed' if failure else 'completed', 'error': failure,
                          'approved_manifest': str(plan_path), 'approved_sha256': args.confirm_manifest_sha256,
                          'retained_verified': retained_verified, 'removed': totals(removed), 'before': context['before_totals'], 'after': totals(after_rows),
                          'accounting': 'File allocation excluding directory blocks and subsequently written retention notes; not a filesystem quota measurement.'}
                write_new(directory / 'execution_summary.json', json_bytes(result))
                for run in runs:
                    count = sum(row['run'] == str(run) for row in removed)
                    note = (f'# Published payload retention\n\nChecked {result["utc"]}. Cleanup status: {result["status"]}.\n'
                            f'Removed {count} published files from this historical run. The original completed-run label no longer implies a full published payload.\n\n'
                            f'Exact manifest: {plan_path}\nSHA-256: {args.confirm_manifest_sha256}\n'
                            f'Per-file outcomes: {directory / "deletions.jsonl"}\n'
                            f'Retained-file verification passed: {retained_verified}. No cache cleanup was performed. '
                            'Recovery is limited to the approved backup/replay policy; Git cannot undo payload deletion.\n')
                    write_new(run / 'PAYLOAD_RETENTION.md', note.encode())
        if failure:
            raise CleanupError(f'Apply stopped; inspect its journal and execution summary: {failure}')
        print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', action='append', required=True, help='Explicit physical direct child of the approved bulk-results root')
    parser.add_argument('--artifacts', default='bam', help='bam (default), context-calls, or bam,context-calls')
    parser.add_argument('--manifest-dir', help='New metadata directory for a dry-run; never overwritten')
    parser.add_argument('--dependency-review', help='Reviewed JSON release evidence, bound into a new dry-run manifest')
    parser.add_argument('--apply', action='store_true', help='Delete only the explicitly approved unchanged manifest')
    parser.add_argument('--approved-manifest')
    parser.add_argument('--confirm-manifest-sha256')
    args = parser.parse_args()
    if args.apply:
        if not args.approved_manifest or not args.confirm_manifest_sha256 or args.manifest_dir or args.dependency_review:
            parser.error('--apply requires --approved-manifest and --confirm-manifest-sha256 only; review is bound inside the approved plan')
    elif not args.manifest_dir or args.approved_manifest or args.confirm_manifest_sha256:
        parser.error('Dry-run requires --manifest-dir; approval options require --apply')
    try:
        apply_plan(args) if args.apply else make_plan(args)
    except (CleanupError, OSError, ValueError, KeyError, TypeError) as exc:
        raise SystemExit(f'ERROR: {exc}') from exc


if __name__ == '__main__':
    main()
