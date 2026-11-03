#!/usr/bin/env python3
"""Audit regenerated inputs against saved cohorts and check native cache recovery."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from common import *
from native import invoke, native_fractions, verify_runtime
from pat_inputs import profile_atlas, read_counts


def validate(root, cache_check=False):
    inventory = read_json(root / 'input_fingerprint.json')['pat_plan']
    names = [r['sample'] for r in inventory['records']]
    jobs = read_json(root / 'inputs/cohort.jobs.json')
    baselines = {b['cohort']: b for b in read_json(WORKFLOW / 'config/baselines.json')}
    differences, masks, counts_verified, native_full = [], [], [], []
    atlas, _ = profile_atlas('full')
    previous_full = Path(read_json(WORKFLOW / 'config/fragment_validation.json')['bundle'])
    for cohort in inventory['settings']['cohorts']:
        baseline = baselines[cohort]
        saved = Path(baseline['epidish'])
        original = read_json(saved / 'provenance/input_manifest.json')
        if original['atlas']['sha256'] != ATLAS_SHA:
            raise ValueError('Historical atlas checksum differs')
        selected = [r['sample'] for r in inventory['records'] if r['cohort'] == cohort]
        if set(selected) != set(original['samples']):
            raise ValueError('Historical validation requires unchanged original cohorts')
        for sample in selected:
            candidates = list((Path(baseline['uxm']) / 'intermediates/l4').glob(f'{sample}.mem.*.homog.gz'))
            if len(candidates) != 1:
                raise ValueError(f'Ambiguous historical count table: {sample}')
            source = next(r for r in original['regional_tables'] if r['sample'] == sample)
            if sha(candidates[0]) != source['sha256']:
                raise ValueError('Historical count checksum changed')
            current = root / 'cohort/counts/fragment_counts' / f'{sample}.uxm.bed.gz'
            np.testing.assert_array_equal(read_counts(current, atlas), read_counts(candidates[0], atlas))
            counts_verified.append(sample)
        for profile in inventory['settings']['profiles']:
            cells = profile_atlas(profile)[1]
            native = native_fractions(root / 'cohort/uxm' / profile, 'uxm', names, cells).loc[selected]
            source = Path(baseline['uxm' if profile == 'full' else 'uxm_hematopoietic']) / 'results/lineage_fractions.csv'
            expected = original['full_uxm_run' if profile == 'full' else 'hematopoietic_uxm_run']['native_result_sha256']
            if sha(source) != expected:
                raise ValueError('Historical UXM fractions changed')
            old = pd.read_csv(source, index_col=0).T.loc[selected, cells]
            delta = (native - old) * 100
            if not np.array_equal(native.values, old.values):
                raise ValueError(f'Native UXM differs from accepted {profile} estimates')
            differences.extend(dict(method='uxm', profile=profile, cohort=cohort, threshold='', sample=s,
                                    cell=c, difference_pp=delta.loc[s, c]) for s in selected for c in cells)
    for job in jobs:
        profile, cells = job['profile'], job['cells']
        historical = Path(baselines[job['cohort']]['epidish']) / f'comparisons/{profile}/coverage-ge{job["threshold"]}'
        manifest = read_json(historical / 'manifest.json')
        files = {f['name']: f['sha256'] for f in manifest['files']}
        for name in ('marker_mask.tsv', 'reference_matrix.tsv', 'sample_regional_values.tsv.gz', 'epidish_rpc_fractions.tsv'):
            if sha(historical / name) != files[name]:
                raise ValueError(f'Historical matrix/fraction checksum changed: {name}')
        fresh_mask = pd.read_csv(Path(job['folder']) / 'marker_mask.tsv', sep='\t')
        old_mask = pd.read_csv(historical / 'marker_mask.tsv', sep='\t')
        if list(fresh_mask.region_id) != list(old_mask.region_id):
            raise ValueError(f'Regenerated marker mask differs: {job["id"]}')
        for name in ('reference_matrix.tsv', 'sample_regional_values.tsv.gz'):
            fresh = pd.read_csv(Path(job['folder']) / name, sep='\t', index_col=0)
            old = pd.read_csv(historical / name, sep='\t', index_col=0).loc[fresh.index, fresh.columns]
            np.testing.assert_allclose(fresh.values, old.values, atol=1e-15, rtol=0)
        coverage_cols = ['minimum_cohort_coverage'] + ['coverage_' + s for s in job['samples']]
        np.testing.assert_array_equal(fresh_mask[coverage_cols].values, old_mask[coverage_cols].values)
        masks.append(dict(id=job['id'], markers=job['markers'], identical=True))
        native = native_fractions(root / 'cohort/epidish' / job['id'], 'epidish', job['samples'], cells)
        old = pd.read_csv(historical / 'epidish_rpc_fractions.tsv', sep='\t', index_col=0).T.loc[job['samples'], cells]
        delta = (native - old) * 100
        differences.extend(dict(method='epidish', profile=profile, cohort=job['cohort'], threshold=job['threshold'],
                                sample=s, cell=c, difference_pp=delta.loc[s, c]) for s in job['samples'] for c in cells)
        if profile == 'full':
            former = native_fractions(previous_full / 'cohort/epidish' / f'{job["cohort"]}_ge{job["threshold"]}',
                                      'epidish', job['samples'], cells)
            error = float((native - former).abs().max().max())
            if error > 1e-8:
                raise ValueError('Regenerated full input changes native fits beyond numeric precision')
            native_full.append(dict(id=job['id'], maximum_absolute_fraction_difference=error))
    frame = pd.DataFrame(differences)
    summary = {f'{m}/{p}': dict(estimates=len(g), mean_absolute_difference_pp=float(g.difference_pp.abs().mean()),
                              maximum_absolute_difference_pp=float(g.difference_pp.abs().max()))
               for (m, p), g in frame.groupby(['method', 'profile'])}
    result = dict(count_tables_identical=counts_verified, marker_masks=masks, comparison=summary,
                  regenerated_vs_saved_native_full=native_full, cache_checks=[])
    if cache_check:
        require_incomplete(root)
        if not os.environ.get('SLURM_JOB_ID'):
            raise ValueError('Submit cache-check through sbatch')
        runtime = verify_runtime(read_json(root / 'input_fingerprint.json')['runtime'])
        targets = [('counts', root / 'inputs/cohort.pats.csv', root / 'inputs/markers.tsv', 'cohort', root / 'cohort/counts'),
                   ('uxm', root / 'inputs/cohort.pats.csv', root / 'inputs/Atlas.U25.l4.hg38.hematopoietic.tsv',
                    'cohort_hematopoietic', root / 'cohort/uxm/hematopoietic')]
        job = next(j for j in jobs if j['profile'] == 'hematopoietic')
        targets.append(('epidish', Path(job['matrix']), Path(job['reference']), 'cohort_' + job['id'], root / 'cohort/epidish' / job['id']))
        frozen = {str(p.relative_to(root)): sha(p) for p in (root / 'cohort').rglob('*')
                  if p.is_file() and p.suffix in ('.csv', '.gz')}
        for method, sheet, reference, tag, dest in targets:
            invoke(runtime, root, method, sheet, reference, tag, resume=True, fragments=True, result_override=dest)
            trace = rows(dest / 'pipeline_info' / f'execution_trace_{tag}.txt')
            if not trace or any(r['status'] != 'CACHED' for r in trace):
                raise ValueError(f'Recovery did not reuse all tasks: {method}')
            result['cache_checks'].append(dict(method=method, cached_tasks=len(trace)))
        for name, expected in frozen.items():
            if sha(root / name) != expected:
                raise ValueError('Recovery changed a completed count/fraction output')
        result['identical_recovery_outputs'] = True
    return result, frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis', type=Path)
    parser.add_argument('--cache-check', action='store_true')
    parser.add_argument('--output', type=Path, required=True, help='New dated validation bundle; source results stay immutable')
    parser.add_argument('--label', default='first')
    args = parser.parse_args()
    root = args.analysis.resolve(strict=True)
    if not (root / 'COHORT_COMPLETE').exists():
        raise ValueError('Native cohort run must complete first')
    if not re.fullmatch(r'[A-Za-z0-9_-]+', args.label):
        raise ValueError('Unsafe validation label')
    if args.cache_check:
        raise ValueError('Completed runs are immutable; cache checks belong to an incomplete execution attempt')
    folder = new_run(args.output)
    target = folder / (args.label + '.json')
    if target.exists():
        raise ValueError('Refusing to overwrite a validation record')
    record = dict(started_utc=utc(), status='running', cache_check=args.cache_check,
                  source_sha256=sha(Path(__file__)), original_execution=read_json(root / 'execution.json'))
    write_json(target, record)
    try:
        result, frame = validate(root, args.cache_check)
        frame.to_csv(folder / (args.label + '.differences.tsv'), sep='\t', index=False)
        record.update(result, status='passed', completed_utc=utc())
    except Exception as exc:
        record.update(status='failed', ended_utc=utc(), error=str(exc))
        raise
    finally:
        write_json(target, record)
        finish(folder, 'complete' if record['status'] == 'passed' else 'failed')


if __name__ == '__main__':
    main()
