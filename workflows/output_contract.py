#!/usr/bin/env python3
"""Output ownership adapter shared by downstream launchers; no scientific defaults."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / 'scripts/lib'))
import output_guard


def identity(workflow, inputs=(), work_root=None):
    workflow = Path(workflow).resolve(strict=True)
    sources = []
    if workflow.is_file():
        sources = [workflow]
    else:
        excluded = {'runs', 'analyses', 'software', 'results', 'outputs', 'exports', 'intermediates',
                    'reference', 'references', '__pycache__', 'pannel_cpg_list', 'logs', 'work', '.nextflow', '.git'}
        for current, directories, filenames in os.walk(workflow, followlinks=False):
            directories[:] = sorted(d for d in directories if d not in excluded and not (Path(current)/d).is_symlink())
            sources.extend(Path(current)/name for name in sorted(filenames)
                if (Path(current)/name).suffix in ('.py', '.R', '.sh', '.config', '.yml', '.yaml')
                and not (Path(current)/name).is_symlink())
    files = [*sources, *(Path(p).resolve(strict=True) for p in inputs)]
    if any(p.stat().st_size > 16 * 1024 * 1024 for p in files):
        raise ValueError('Ownership hashes accept small source/configuration manifests only; use recorded bulk identities')
    revision = output_guard.require_clean_source(PROJECT)
    source_revision = output_guard.require_clean_source(workflow)
    return {'workflow': str(workflow), 'project_revision': revision, 'source_revision': source_revision,
            'submitted_inputs': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'work_root': str(Path(work_root).resolve()) if work_root else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=('reserve', 'resume', 'check', 'complete', 'fail', 'publish'))
    parser.add_argument('output', type=Path)
    parser.add_argument('--workflow', type=Path)
    parser.add_argument('--input', type=Path, action='append', default=[])
    parser.add_argument('--work-root', type=Path)
    parser.add_argument('--source', type=Path)
    args = parser.parse_args()
    if args.operation in ('reserve', 'resume'):
        if not args.workflow: parser.error('--workflow required')
        stamp = identity(args.workflow, args.input, args.work_root)
        if args.operation == 'reserve':
            output_guard.reserve_path(args.output, stamp)
            output_guard.update(args.output, status='running')
        else: output_guard.resume(args.output, stamp)
    elif args.operation == 'check':
        revision = output_guard.require_clean_source(PROJECT)
        record = output_guard.read_record(args.output)
        if record['identity'].get('project_revision') != revision:
            raise ValueError('Execution requires the original project revision')
        if record.get('status') in ('complete', 'completed'):
            raise ValueError('Completed output is immutable')
    elif args.operation in ('complete', 'fail'):
        output_guard.finish(args.output, args.operation == 'complete')
    elif args.operation == 'publish':
        if not args.source: parser.error('--source required')
        output_guard.publish_contents(args.source, args.output)


if __name__ == '__main__':
    main()
