#!/usr/bin/env python3
"""Output ownership adapter for relocated experimental entrypoints.

No resume is implied by this adapter. New bundles and internal sections are
exclusive creations; historical runs continue through their frozen launchers.
"""
import argparse
import functools
import hashlib
import inspect
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
from datetime import datetime, timezone

import output_guard as guard

_reserved = []


def utc():
    return datetime.now(timezone.utc).isoformat()


def revision(path):
    result = subprocess.run(['git', '-C', str(path), 'rev-parse', 'HEAD'],
                            capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else None


def require_current_sources(source=None):
    """Require committed execution helpers and, for new work, its source repo."""
    project = Path(__file__).resolve().parents[2]
    project_revision = guard.require_clean_source(project)
    experiment_revision = guard.require_clean_source(Path(source).resolve().parent) if source else None
    return project_revision, experiment_revision


def source_identity(source, inputs=()):
    source = Path(source).resolve()
    project_revision, experiment_revision = require_current_sources(source)
    identities = {}
    for item in inputs:
        path = Path(item)
        if not path.exists():
            identities[str(path)] = {'missing': True}
            continue
        st = path.stat()
        value = dict(path=str(path.resolve()), size=st.st_size, mtime_ns=st.st_mtime_ns)
        if path.is_file() and st.st_size <= 1024 * 1024:
            value['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        identities[str(path)] = value
    return dict(source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                experiment_revision=experiment_revision, project_revision=project_revision,
                guard_sha256=hashlib.sha256(Path(guard.__file__).read_bytes()).hexdigest(),
                adapter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                inputs=identities, resume_supported=False)


def _source_and_inputs():
    # Capture the invoking script and explicit argparse Path inputs without
    # scanning or hashing bulk data. Scientific manifests remain authoritative.
    first = None
    for frame in inspect.stack()[2:]:
        if Path(frame.filename).resolve() == Path(__file__).resolve():
            continue
        first = first or frame.filename
        args = frame.frame.f_locals.get('args', frame.frame.f_locals.get('a'))
        if not isinstance(args, argparse.Namespace):
            continue
        values = vars(args).items() if isinstance(args, argparse.Namespace) else ()
        inputs = []
        for key, value in values:
            if any(word in key for word in ('output', 'out_dir', 'run_dir', 'run_root', 'tmp', 'cache')):
                continue
            for item in value if isinstance(value, (list, tuple)) else [value]:
                if isinstance(item, Path) or (isinstance(item, str) and Path(item).is_file()):
                    inputs.append(item)
        return frame.filename, inputs
    if first:
        return first, []
    raise RuntimeError('Cannot identify experiment source')


def _canonical(path):
    if '..' in Path(path).parts:
        raise ValueError('Output traversal is forbidden')
    path = Path(os.path.abspath(path))
    if '..' in Path(path).parts or path.parent.resolve() != path.parent:
        raise ValueError('Output parent must be canonical, without symlink aliases')
    if os.path.lexists(path):
        raise ValueError(f'Output already exists; choose a new directory: {path}')
    _writable_ancestors(path.parent)
    return path


def _writable_ancestors(path):
    """Check the actual destination tree before making even intermediate dirs."""
    path = Path(path)
    while not os.path.lexists(path):
        path = path.parent
    guard.canonical_parent(path)


def reserve_bundle(output, run_root=None, run_name=None):
    output = _canonical(output)
    source, inputs = _source_and_inputs()
    identity = source_identity(source, inputs)
    if run_root is None and run_name is None:
        reserved = guard.reserve(output.parent, output.name, identity)
        guard.update(reserved, status='running')
    else:
        # Internal sections may have simple names. They must be below either a
        # real incomplete owned run or the launcher's exact staging directory.
        if run_root is not None:
            root = Path(run_root)
            if root.is_symlink() or root.resolve() != root.absolute():
                raise ValueError('Section run root must be canonical')
            record = guard.read_record(root)
            if record.get('status') in {'completed', 'complete'}:
                raise ValueError('Completed run cannot receive a new section')
            output.relative_to(root.resolve(strict=True))
            owner = root
        else:
            from datetime import datetime
            if not re.fullmatch(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-[A-Za-z0-9][A-Za-z0-9._-]*', run_name):
                raise ValueError('Invalid staging run name')
            datetime.strptime(run_name[:16], '%Y-%m-%d_%H-%M')
            matches = [p for p in output.parents if re.fullmatch(
                r'\.' + re.escape(run_name) + r'\.running-[A-Za-z0-9._-]+', p.name)]
            if not matches:
                raise ValueError('Staging output does not match explicit run name')
            owner = os.environ.get('EXPERIMENT_OWNER_OUTPUT')
            if not owner or Path(owner).name != run_name:
                raise ValueError('Staging requires its reserved final bundle')
            record = guard.read_record(owner)
            if record.get('status') in {'completed', 'complete'}:
                raise ValueError('Completed run cannot receive a new section')
        output.parent.mkdir(parents=True, exist_ok=True)
        os.mkdir(output)
        # Sections are independently exclusive; the final run remains owned by
        # its shell launcher. Their receipt intentionally has no resume route.
        with (output / '.run-ownership.json').open('x') as handle:
            json.dump(dict(schema=1, run_id=output.name, output=str(output), identity=identity,
                           started_utc=utc(), status='running', section=True,
                           owner_output=str(Path(owner).resolve()),
                           session_id=None, attempts=[]), handle, indent=2)
        reserved = output
    _reserved.append(reserved)
    return reserved


def reserve_resource(output, source=None):
    """References are immutable named resources, exempt from dated run names."""
    output = _canonical(output)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', output.name):
        raise ValueError('Invalid resource name')
    source, inputs = (source, []) if source else _source_and_inputs()
    identity = source_identity(source, inputs)
    guard.canonical_parent(output.parent)
    os.mkdir(output)
    with (output / '.run-ownership.json').open('x') as handle:
        json.dump(dict(schema=1, run_id=output.name, output=str(output),
                       identity=identity, status='running',
                       started_utc=utc(), session_id=None, attempts=[], resource=True), handle, indent=2)
    _reserved.append(output)
    return output


def freeze_prepared_run(output):
    """Bind a multi-stage experiment to its exact immutable submission files."""
    output = Path(output)
    guard.update(output,
                 snapshot_manifest_sha256=hashlib.sha256((output / 'source_hashes.json').read_bytes()).hexdigest(),
                 input_manifest_sha256=hashlib.sha256((output / 'inputs.json').read_bytes()).hexdigest())


def require_owned_run(output, source=None):
    require_current_sources()
    record = guard.read_record(output)
    if record['status'] not in {'running', 'failed'}:
        raise ValueError('Run is complete or was never initialized; use a new run name')
    for field, name in [('snapshot_manifest_sha256', 'source_hashes.json'), ('input_manifest_sha256', 'inputs.json')]:
        actual = hashlib.sha256((Path(output) / name).read_bytes()).hexdigest()
        if record.get(field) != actual:
            raise ValueError('Frozen run identity differs or is missing: ' + name)
    for name, expected in json.loads((Path(output) / 'source_hashes.json').read_text()).items():
        actual = hashlib.sha256((Path(output) / 'submitted' / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Submitted source differs: ' + name)
    if record['identity'].get('guard_sha256') != hashlib.sha256(Path(guard.__file__).read_bytes()).hexdigest():
        raise ValueError('Output guard differs from prepared run')
    if record['identity'].get('adapter_sha256') != hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
        raise ValueError('Experiment adapter differs from prepared run')
    if source is not None:
        source = Path(source).resolve()
        root = Path(output).resolve()
        try:
            relative = source.relative_to(root)
        except ValueError:
            # A current recovery/launch helper must itself be committed.
            require_current_sources(source)
        else:
            expected = None
            if relative.parts[:1] == ('submitted',):
                expected = json.loads((root / 'source_hashes.json').read_text()).get(
                    relative.relative_to('submitted').as_posix())
            elif relative.parts[:1] == ('recovery',) and len(relative.parts) >= 3:
                attempt = root / Path(*relative.parts[:2])
                manifest = json.loads((attempt / 'manifest.json').read_text())
                if Path(manifest['run_dir']).resolve() != root:
                    raise ValueError('Recovery manifest belongs to another run')
                expected = manifest['recovery_source_hashes'].get(source.relative_to(attempt).as_posix())
            if expected != hashlib.sha256(source.read_bytes()).hexdigest():
                raise ValueError('Executing source is not a verified submitted snapshot')
    return record


def require_incomplete_parent(output):
    """Internal tools may write only beneath a canonical incomplete bundle."""
    output = Path(os.path.abspath(output))
    if output.resolve() != output:
        raise ValueError('Internal output cannot use symlink aliases')
    _writable_ancestors(output)
    for ancestor in [output, *output.parents]:
        if (ancestor / '.run-ownership.json').exists():
            record = guard.read_record(ancestor)
            if record['status'] != 'running':
                raise ValueError('Internal writes require a running owned bundle')
            return ancestor
    raise ValueError('Internal output has no owned run ancestor')


def prepare_file_outputs(paths):
    """Direct reports reserve a bundle; launcher reports stay inside its owner."""
    source, _ = _source_and_inputs()
    require_current_sources(source)
    paths = [Path(p).absolute() for p in paths if p is not None]
    for path in paths:
        if os.path.lexists(path):
            raise ValueError('Output file already exists: ' + str(path))
        if path.parent.resolve() != path.parent:
            raise ValueError('Output file parent cannot use an alias')
        _writable_ancestors(path.parent)
    owner = os.environ.get('EXPERIMENT_OWNER_OUTPUT')
    if owner:
        owner = require_incomplete_parent(owner)
        for path in paths:
            path.relative_to(owner)
    else:
        parents = {p.parent for p in paths}
        if len(parents) != 1:
            raise ValueError('Direct report files must share a new dated bundle')
        reserve_bundle(parents.pop())


def finish_bundle(output, success):
    record_path = Path(output) / '.run-ownership.json'
    if not record_path.exists():
        return
    guard.finish(output, success=success)


def owned_execution(function):
    """Finalize only bundles reserved by this call; failures retain evidence."""
    @functools.wraps(function)
    def wrapped(*args, **kwargs):
        offset = len(_reserved)
        previous = signal.getsignal(signal.SIGTERM)
        def interrupted(signum, frame):
            raise SystemExit(128 + signum)
        signal.signal(signal.SIGTERM, interrupted)
        try:
            result = function(*args, **kwargs)
        except BaseException:
            for output in reversed(_reserved[offset:]):
                if guard.read_record(output)['status'] != 'completed':
                    finish_bundle(output, False)
            raise
        else:
            for output in reversed(_reserved[offset:]):
                if guard.read_record(output)['status'] != 'completed':
                    finish_bundle(output, result in (None, 0))
            return result
        finally:
            signal.signal(signal.SIGTERM, previous)
    return wrapped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['reserve', 'reserve-resource', 'publish', 'finish', 'link', 'check-parent', 'claim-file'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--input', type=Path, action='append', default=[])
    parser.add_argument('--staging', type=Path)
    parser.add_argument('--link', type=Path)
    parser.add_argument('--success', action='store_true')
    args = parser.parse_args()
    require_current_sources()
    if args.action in {'reserve', 'reserve-resource'}:
        if args.source is None:
            parser.error('--source is required')
        if args.action == 'reserve-resource':
            reserve_resource(args.output, args.source)
        else:
            guard.reserve(args.output.parent, args.output.name, source_identity(args.source, args.input))
            guard.update(args.output, status='running')
    elif args.action == 'publish':
        if args.staging is None:
            parser.error('--staging is required')
        destination = guard.read_record(args.output)
        if destination['status'] != 'running':
            raise ValueError('Publication requires a running destination owner')
        staging = Path(os.path.abspath(args.staging))
        if staging.resolve(strict=True) != staging or not staging.is_dir():
            raise ValueError('Publication staging must be a canonical directory')
        if not any(re.fullmatch(r'\.' + re.escape(Path(args.output).name) +
                               r'\.running-[A-Za-z0-9._-]+', p.name)
                   for p in [staging, *staging.parents]):
            raise ValueError('Publication source is not the matching run staging directory')
        _writable_ancestors(staging.parent)
        record = args.staging / '.run-ownership.json'
        archived_record = args.staging / 'staging-ownership.json'
        section = None
        if os.path.lexists(record):
            section = guard.read_record(staging)
        elif os.path.lexists(archived_record):
            if archived_record.is_symlink():
                raise ValueError('Staging ownership receipt cannot be a symlink')
            section = json.loads(archived_record.read_text())
        if section is not None:
            if not section.get('section') or section.get('owner_output') != destination['output']:
                raise ValueError('Publication source is not this owner\'s internal section')
            if section['status'] != 'completed':
                raise ValueError('Internal section has not completed')
        for child in staging.iterdir():
            name = 'staging-ownership.json' if child.name == '.run-ownership.json' else child.name
            if os.path.lexists(args.output / name):
                raise FileExistsError('Publication conflict: ' + name)
        if os.path.lexists(record):
            # Preserve the internal-section receipt under its own name. The
            # final bundle retains the launcher's original ownership record.
            destination = args.staging / 'staging-ownership.json'
            with destination.open('x') as handle:
                handle.write(record.read_text())
            record.unlink()
        guard.publish_contents(args.staging, args.output)
    elif args.action == 'finish':
        guard.finish(args.output, success=args.success)
    elif args.action == 'check-parent':
        require_incomplete_parent(args.output)
    elif args.action == 'claim-file':
        require_incomplete_parent(args.output.parent)
        with _canonical(args.output).open('x'):
            pass
    else:
        if args.link is None:
            parser.error('--link is required')
        guard.publish_link(args.link, args.output)


if __name__ == '__main__':
    main()
