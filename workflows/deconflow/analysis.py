#!/usr/bin/env python3
"""Run native UXM from an explicit PAT plan; patient CelFiE remains reference-gated."""
import argparse
import fcntl
import os
from common import *
from native import verify_runtime, invoke, native_fractions
from pat_inputs import read_plan, prepare_pat_inputs, prepared_hashes, profile_atlas


def selected_samples(data, phase):
    if len(data) != 12 or len({r['sample'] for r in data}) != 12:
        raise ValueError('Expected twelve unique sample identifiers')
    patients = list(dict.fromkeys(r['patient'] for r in data))
    if len(patients) != 4:
        raise ValueError('Expected four patients')
    for patient in patients:
        group = [r for r in data if r['patient'] == patient]
        if len(group) != 3 or {r['compartment'] for r in group} != {'BM', 'PB', 'cfDNA'}:
            raise ValueError('Incomplete or duplicated patient triplet')
    return data if phase == 'cohort' else [r for r in data if r['patient'] == patients[0]]


def code_fingerprint():
    return {str(p.relative_to(WORKFLOW)): sha(p) for p in
            list(WORKFLOW.glob('*.py')) + list(WORKFLOW.glob('*.sh')) + list(WORKFLOW.glob('*.nf')) +
            list((WORKFLOW / 'config').glob('*.config')) + list((WORKFLOW / 'config').glob('params.*.json'))}


def run_pat_plan(args, runtime):
    if args.phase != 'cohort':
        raise ValueError('PAT plans use --phase cohort; select the intended samples explicitly in the plan')
    inventory = read_plan(args.plan)
    current = dict(runtime=str(runtime), runtime_lock=sha(runtime / 'runtime.lock.json'),
                   pat_plan=inventory, code=code_fingerprint())
    root = args.output.absolute()
    if args.resume:
        require_incomplete(root)
        if not root.is_dir() or read_json(root / 'input_fingerprint.json') != current:
            raise ValueError('Resume requires unchanged PAT plan, input identities, atlas, runtime and code')
        if (root / 'COHORT_COMPLETE').exists():
            raise ValueError('Completed cohort cannot be overwritten')
    else:
        root = new_run(root)
    with (root / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        require_incomplete(root)
        if not args.resume:
            write_json(root / 'input_fingerprint.json', current)
            prepare_pat_inputs(root, inventory)
            write_json(root / 'prepared.sha256.json', prepared_hashes(root))
        attempt = dict(job=os.environ['SLURM_JOB_ID'], phase='cohort', started_utc=utc(), resume=args.resume)
        try:
            for name, expected in read_json(root / 'prepared.sha256.json').items():
                if sha(root / name) != expected:
                    raise ValueError(f'Prepared UXM input changed: {name}')
            wanted = [r['sample'] for r in inventory['records']]
            for profile in inventory['settings']['profiles']:
                result = root / 'cohort/uxm' / profile
                invoke(runtime, root, 'uxm', root / 'inputs/cohort.pats.csv',
                       root / 'inputs' / f'Atlas.U25.l4.hg38.{profile}.tsv', f'cohort_{profile}',
                       args.resume, fragments=True, result_override=result)
                native_fractions(result, 'uxm', wanted, profile_atlas(profile)[1])
            (root / 'COHORT_COMPLETE').write_text(utc() + '\n')
            attempt.update(status='complete', completed_utc=utc())
            finish(root, phase='cohort', profiles=inventory['settings']['profiles'], samples=wanted,
                   methods=['uxm'], input_representation='indexed hg38 PATs; original U25/l4 atlas',
                   fitting='unchanged DecoNFlow UXM module')
        except Exception as exc:
            attempt.update(status='failed', ended_utc=utc(), error=str(exc))
            finish(root, 'failed', error=str(exc))
            raise
        finally:
            with (root / 'attempts.jsonl').open('a') as handle:
                handle.write(json.dumps(attempt) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--validation', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True,
                        help='Explicit PAT samplesheet and profiles; historical cohort thresholds are retained but unused by UXM')
    parser.add_argument('--phase', choices=['triplet', 'cohort'], default='triplet')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    if args.phase != 'cohort':
        parser.error('Use --phase cohort with an explicit PAT plan; the historical UXM/EpiDISH route is retired')
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Use run_deconflow.sh through sbatch')
    runtime = verify_runtime(args.runtime)
    gate = args.validation / 'VALIDATED.json'
    if not gate.is_file() or read_json(gate)['runtime_lock'] != sha(runtime / 'runtime.lock.json'):
        raise ValueError('Native offline/runtime/recovery fixtures have not passed for this runtime')
    run_pat_plan(args, runtime)


if __name__ == '__main__':
    main()
