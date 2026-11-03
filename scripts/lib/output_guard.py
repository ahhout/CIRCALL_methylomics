#!/usr/bin/env python3
"""Exclusive result-bundle ownership. Only caller-owned incomplete runs may change."""
import argparse
import ctypes
import errno
import json
import os
import re
import stat
import tempfile
from datetime import datetime, timezone
from pathlib import Path

RECORD = '.run-ownership.json'


def utc_now():
    return datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')


def canonical_parent(parent):
    parent = Path(os.path.abspath(os.fspath(parent)))
    resolved = parent.resolve(strict=True)
    if parent != resolved or not resolved.is_dir():
        raise ValueError('Results parent must be an existing canonical directory, not a symlink alias')
    for ancestor in (resolved, *resolved.parents):
        marker = ancestor / RECORD
        if os.path.lexists(marker):
            if marker.is_symlink() or not marker.is_file():
                raise ValueError('Invalid ancestor ownership marker')
            with marker.open() as handle:
                owner = json.load(handle)
            if owner.get('status') == 'completed':
                raise ValueError('Cannot write inside a completed result bundle')
        if (ancestor / 'COMPLETED').is_file():
            raise ValueError('Cannot write inside a completed historical result bundle')
    return resolved


def validate_name(name):
    if not isinstance(name, str) or not re.fullmatch(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-[A-Za-z0-9][A-Za-z0-9_.-]*', name):
        raise ValueError('Run name must be a direct child named YYYY-MM-DD_HH-MM-description (UTC)')
    datetime.strptime(name[:16], '%Y-%m-%d_%H-%M')
    return name


def _output(output):
    raw = os.fspath(output)
    if '..' in Path(raw).parts:
        raise ValueError('Output traversal is forbidden')
    p = Path(os.path.abspath(raw))
    parent = canonical_parent(p.parent)
    if p.is_symlink() or not p.is_dir() or p.resolve(strict=True) != p:
        raise ValueError('Output must be an owned regular directory, not an alias')
    if p.parent != parent:
        raise ValueError('Output escapes its parent')
    return p


def _write(output, record):
    fd, temp = tempfile.mkstemp(prefix='.ownership-', dir=output)
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(record, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        target = output / RECORD
        if target.is_symlink():
            raise ValueError('Ownership record cannot be a symlink')
        os.replace(temp, target)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def reserve(parent, run_name, identity):
    """Reserve final leaf atomically; never adopt even an empty existing directory."""
    validate_name(run_name)
    if not isinstance(identity, dict) or not identity:
        raise ValueError('A nonempty identity object is required')
    # Validate serializability before reserving any path.
    json.dumps(identity, allow_nan=False)
    if '..' in Path(os.fspath(parent)).parts:
        raise ValueError('Parent traversal is forbidden')
    p = canonical_parent(parent) / run_name
    os.mkdir(p)  # EEXIST includes files and dangling symlinks; do not remove on failure.
    record = dict(schema=1, run_id=run_name, output=str(p), identity=identity,
                  status='reserved', started_utc=utc_now(), session_id=None, attempts=[])
    _write(p, record)  # A failed initialization deliberately leaves the reserved directory.
    return p


def reserve_path(output, identity):
    p = Path(output)
    if '..' in p.parts:
        raise ValueError('Output traversal is forbidden')
    return reserve(p.parent, p.name, identity)


def read_record(output):
    p = _output(output)
    fd = os.open(p / RECORD, os.O_RDONLY | os.O_NOFOLLOW)
    with os.fdopen(fd) as handle:
        if not stat.S_ISREG(os.fstat(handle.fileno()).st_mode):
            raise ValueError('Ownership record is not a regular file')
        record = json.load(handle)
    required = {'schema', 'run_id', 'output', 'identity', 'status', 'started_utc', 'session_id', 'attempts'}
    if not required.issubset(record) or record['schema'] != 1 or record['output'] != str(p) or record['run_id'] != p.name:
        raise ValueError('Invalid ownership record or output identity')
    if record['status'] not in {'reserved', 'running', 'failed', 'completed'}:
        raise ValueError('Invalid ownership status')
    if not isinstance(record['identity'], dict) or not isinstance(record['attempts'], list):
        raise ValueError('Invalid ownership identity or attempts')
    return record


def resume(output, identity, session_id=None):
    record = read_record(output)
    if record['status'] not in {'running', 'failed'}:
        raise ValueError('Only an explicitly owned incomplete started run may resume')
    if record['identity'] != identity:
        raise ValueError('Resume identity differs from original submitted inputs/source')
    if session_id is not None and record['session_id'] != session_id:
        raise ValueError('Resume session identity differs')
    return record


def update(output, **fields):
    p = _output(output)
    record = read_record(p)
    if record['status'] == 'completed':
        raise ValueError('Completed outputs are immutable')
    if {'schema', 'run_id', 'output', 'identity', 'started_utc'}.intersection(fields):
        raise ValueError('Original ownership identity cannot change')
    if 'session_id' in fields and record['session_id'] is not None and fields['session_id'] != record['session_id']:
        raise ValueError('Original session identity cannot change')
    if 'status' in fields and fields['status'] not in {'reserved', 'running', 'failed', 'completed'}:
        raise ValueError('Invalid status')
    record.update(fields)
    _write(p, record)
    return record


def finish(output, success, **fields):
    fields.update(status='completed' if success else 'failed')
    fields['completed_utc' if success else 'ended_utc'] = utc_now()
    return update(output, **fields)


def publish_link(link, output):
    """Create a navigation link exclusively, accepting only an identical existing link."""
    p = _output(output)
    read_record(p)
    link = Path(link)
    if not link.parent.is_dir():
        raise ValueError('Navigation parent must already exist')
    canonical_parent(link.parent)
    try:
        os.symlink(str(p), link)
    except FileExistsError:
        if not link.is_symlink() or link.resolve(strict=False) != p:
            raise ValueError('Navigation path already belongs to another output')


def _rename_noreplace(source, destination):
    libc = ctypes.CDLL(None, use_errno=True)
    rename = getattr(libc, 'renameat2', None)
    if rename is None:
        raise OSError(errno.ENOSYS, 'renameat2 unavailable')
    rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    rename.restype = ctypes.c_int
    if rename(-100, os.fsencode(source), -100, os.fsencode(destination), 1):
        code = ctypes.get_errno()
        raise OSError(code, os.strerror(code), str(source))


def _link_transfer(source, destination):
    """BeeGFS fallback: exclusive metadata operations, no payload copying/replacement."""
    info = source.lstat()
    if stat.S_ISLNK(info.st_mode):
        os.symlink(os.readlink(source), destination)
        source.unlink()
    elif stat.S_ISREG(info.st_mode):
        os.link(source, destination, follow_symlinks=False)
        source.unlink()
    elif stat.S_ISDIR(info.st_mode):
        os.mkdir(destination, stat.S_IMODE(info.st_mode))
        for child in sorted(source.iterdir()):
            _link_transfer(child, destination / child.name)
        os.utime(destination, ns=(info.st_atime_ns, info.st_mtime_ns))
        # Retain empty source directories as recovery evidence.
    else:
        raise ValueError('Unsupported staged file type')


def publish_contents(source, output):
    """Publish staged children exclusively; retain evidence on any partial failure."""
    p = _output(output)
    record = read_record(p)
    if record['status'] == 'completed':
        raise ValueError('Completed outputs are immutable')
    raw_source = Path(source)
    source = Path(os.path.abspath(raw_source))
    if ('..' in raw_source.parts or source.is_symlink() or not source.is_dir()
            or source.resolve(strict=True) != source or source == p):
        raise ValueError('Invalid staging directory')
    # Publication removes source entries as well as writing destination entries.
    # Protect completed source ancestors before moving even the first child.
    canonical_parent(source.parent)
    if os.path.lexists(source / RECORD):
        if read_record(source)['status'] == 'completed':
            raise ValueError('Cannot publish from a completed result bundle')
    receipt = source / 'staging-ownership.json'
    if os.path.lexists(receipt):
        if receipt.is_symlink() or not receipt.is_file():
            raise ValueError('Invalid staging ownership receipt')
        section = json.loads(receipt.read_text())
        if (section.get('section') is not True or section.get('owner_output') != str(p)
                or not any(path.name.startswith('.' + p.name + '.running-')
                           for path in (source, *source.parents))):
            raise ValueError('Staging receipt does not belong to this output')
    elif os.path.lexists(source / 'COMPLETED'):
        raise ValueError('Cannot publish from a completed historical result bundle')
    entries = sorted(source.iterdir())
    for child in entries:
        if child.name == RECORD or os.path.lexists(p / child.name):
            raise FileExistsError('Staged child conflicts with reserved destination: ' + child.name)
    for child in entries:
        try:
            _rename_noreplace(child, p / child.name)
        except OSError as exc:
            if exc.errno not in {errno.ENOSYS, errno.EINVAL, errno.EOPNOTSUPP}:
                raise
            _link_transfer(child, p / child.name)
    # Keep the now-empty staging directory as recovery evidence.



def require_clean_source(path):
    """Require an existing committed repository with no tracked or untracked edits."""
    import subprocess
    source = Path(path).resolve(strict=True)
    directory = source if source.is_dir() else source.parent
    try:
        root = subprocess.check_output(
            ['git', '-C', str(directory), 'rev-parse', '--show-toplevel'],
            text=True, stderr=subprocess.PIPE).strip()
        revision = subprocess.check_output(
            ['git', '-C', root, 'rev-parse', '--verify', 'HEAD'],
            text=True, stderr=subprocess.PIPE).strip()
        status = subprocess.check_output(
            ['git', '-C', root, 'status', '--porcelain=v1', '--untracked-files=normal'],
            text=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as exc:
        raise ValueError('Execution requires committed source in a Git repository') from exc
    if not re.fullmatch(r'[0-9a-f]{40,64}', revision):
        raise ValueError('Execution requires a valid committed source revision')
    if status:
        raise ValueError('Execution requires clean committed source: ' + root)
    return revision


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('reserve'); p.add_argument('--parent', required=True); p.add_argument('--name', required=True); p.add_argument('--identity', required=True)
    p = sub.add_parser('resume'); p.add_argument('--output', required=True); p.add_argument('--identity', required=True); p.add_argument('--session-id')
    p = sub.add_parser('finish'); p.add_argument('--output', required=True); p.add_argument('--success', action='store_true')
    p = sub.add_parser('update'); p.add_argument('--output', required=True); p.add_argument('--fields', required=True)
    p = sub.add_parser('link'); p.add_argument('--link', required=True); p.add_argument('--output', required=True)
    p = sub.add_parser('publish'); p.add_argument('--source', required=True); p.add_argument('--output', required=True)
    p = sub.add_parser('record'); p.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.command == 'reserve': print(reserve(args.parent, args.name, json.loads(args.identity)))
    elif args.command == 'resume': print(json.dumps(resume(args.output, json.loads(args.identity), args.session_id)))
    elif args.command == 'finish': print(json.dumps(finish(args.output, args.success)))
    elif args.command == 'update': print(json.dumps(update(args.output, **json.loads(args.fields))))
    elif args.command == 'link': publish_link(args.link, args.output)
    elif args.command == 'publish': publish_contents(args.source, args.output)
    else: print(json.dumps(read_record(args.output)))


if __name__ == '__main__':
    main()
