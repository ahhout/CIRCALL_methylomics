"""Validate saved fragment matrices and adapt file layout, without refitting."""
import numpy as np
import pandas as pd
from common import *


def write_matrices(dest, tag, mask, ref, samples):
    """Adapt row labels to the unchanged native wrapper's coordinate interface."""
    ref, samples = ref.copy(), samples.copy()
    ref.index = mask['name']
    samples.index = mask['name']
    matrix = pd.concat([mask[['chr', 'start', 'end']].set_axis(samples.index), samples], axis=1)
    ref_path, sample_path = dest / f'{tag}.ref.csv', dest / f'{tag}.samples.csv'
    ref.to_csv(ref_path)
    matrix.to_csv(sample_path)
    return dict(reference=str(ref_path), matrix=str(sample_path))


def saved_inputs():
    atlas, cells = atlas_rows()
    canonical = pd.DataFrame(atlas).set_index('name')[cells].apply(pd.to_numeric, errors='coerce')
    pats, jobs, sources = {}, [], {}
    for baseline in read_json(WORKFLOW / 'config/baselines.json'):
        epi, uxm = Path(baseline['epidish']), Path(baseline['uxm'])
        manifest_path = epi / 'provenance/input_manifest.json'
        manifest = read_json(manifest_path)
        sources[str(manifest_path)] = sha(manifest_path)
        if manifest['atlas']['sha256'] != ATLAS_SHA:
            raise ValueError('Saved fragment atlas differs from the canonical UXM atlas')
        pat_manifest = uxm / 'pat_manifest.tsv'
        sources[str(pat_manifest)] = sha(pat_manifest)
        for row in rows(pat_manifest):
            if row['sample'] in pats:
                raise ValueError('Overlapping historical PAT cohorts')
            pats[row['sample']] = {key: identity(uxm / row[key]) for key in ('pat', 'csi')}
        if sha(uxm / 'results/lineage_fractions.csv') != manifest['full_uxm_run']['native_result_sha256']:
            raise ValueError('Saved UXM fractions changed')
        sources[str(uxm / 'results/lineage_fractions.csv')] = sha(uxm / 'results/lineage_fractions.csv')
        for threshold in manifest['coverage_thresholds']:
            folder = epi / f'comparisons/full/coverage-ge{threshold}'
            record = read_json(folder / 'manifest.json')
            sources[str(folder / 'manifest.json')] = sha(folder / 'manifest.json')
            required = {'marker_mask.tsv', 'reference_matrix.tsv', 'sample_regional_values.tsv.gz',
                        'epidish_rpc_fractions.tsv', 'epidish_warnings.tsv'}
            files = {r['name']: r['sha256'] for r in record['files']}
            if not required <= files.keys():
                raise ValueError('Incomplete saved fragment manifest')
            for name in sorted(required):
                path = folder / name
                if sha(path) != files[name]:
                    raise ValueError(f'Saved fragment checksum changed: {path}')
                sources[str(path)] = files[name]
            mask = pd.read_csv(folder / 'marker_mask.tsv', sep='\t')
            ref = pd.read_csv(folder / 'reference_matrix.tsv', sep='\t', index_col=0)
            sample = pd.read_csv(folder / 'sample_regional_values.tsv.gz', sep='\t', index_col=0)
            if (record['profile'] != 'full' or record['coverage_threshold'] != threshold or
                    record['samples'] != manifest['samples'] or list(ref.columns) != cells or
                    list(sample.columns) != record['samples'] or len(ref) != record['marker_count'] or
                    list(ref.index) != list(mask.region_id) or list(sample.index) != list(ref.index) or
                    mask['name'].duplicated().any()):
                raise ValueError('Saved fragment labels, dimensions or marker order changed')
            if not np.allclose(ref.values, canonical.loc[mask['name']].values, atol=1e-15, rtol=0):
                raise ValueError('Saved reference values differ from canonical atlas scores')
            for frame in (ref, sample):
                if not np.isfinite(frame.values).all() or ((frame.values < 0) | (frame.values > 1)).any():
                    raise ValueError('Saved fragment matrix contains missing or invalid values')
            if not (mask.minimum_cohort_coverage >= threshold).all():
                raise ValueError('Saved cohort coverage mask changed')
            jobs.append(dict(id=f"{baseline['cohort']}_ge{threshold}", cohort=baseline['cohort'],
                             threshold=threshold, folder=str(folder), samples=record['samples'],
                             markers=len(ref), historical_maxit=manifest['epidish_maxit'],
                             historical_epidish_version=manifest['epidish_version']))
    expected = {r['sample'] for r in rows(WORKFLOW / 'config/samples.tsv')}
    if set(pats) != expected or set(s for job in jobs for s in job['samples']) != expected:
        raise ValueError('Historical fragment inputs do not cover exactly the twelve samples')
    return dict(atlas_sha256=ATLAS_SHA, sources=sources, pats=pats, jobs=jobs)


def prepare_fragments(root, records, inventory):
    dest = root / 'inputs'
    dest.mkdir()
    for sample, artifacts in inventory['pats'].items():
        for key, suffix in [('pat', '.pat.gz'), ('csi', '.pat.gz.csi')]:
            (dest / (sample + suffix)).symlink_to(artifacts[key]['path'])
    from analysis import selected_samples
    for phase in ('triplet', 'cohort'):
        selected = [r['sample'] for r in selected_samples(records, phase)]
        table(dest / f'{phase}.pats.csv', [dict(sample=s, pat=str(dest / (s + '.pat.gz')),
              csi=str(dest / (s + '.pat.gz.csi'))) for s in selected], ['sample', 'pat', 'csi'], ',')
        jobs = []
        for job in inventory['jobs']:
            wanted = [s for s in selected if s in job['samples']]
            if not wanted:
                continue
            folder = Path(job['folder'])
            mask = pd.read_csv(folder / 'marker_mask.tsv', sep='\t')
            ref = pd.read_csv(folder / 'reference_matrix.tsv', sep='\t', index_col=0)
            samples = pd.read_csv(folder / 'sample_regional_values.tsv.gz', sep='\t', index_col=0)[wanted]
            paths = write_matrices(dest, f'{phase}_{job["id"]}', mask, ref, samples)
            jobs.append(dict(job, samples=wanted, **paths))
        write_json(dest / f'{phase}.jobs.json', jobs)
    shutil.copy2(WORKFLOW / 'config/samples.tsv', root / 'samples.tsv')
    shutil.copy2(WORKFLOW / 'config/baselines.json', root / 'baselines.json')
    write_json(root / 'prepared.sha256.json', {str(p.relative_to(root)): sha(p) for p in dest.iterdir() if p.is_file() and not p.is_symlink()})


def compare_fragments(root, phase, cells):
    from native import native_fractions, checked
    differences = []
    baselines = read_json(root / 'baselines.json')
    selected = [r['sample'] for r in __import__('analysis').selected_samples(rows(root / 'samples.tsv'), phase)]
    native = native_fractions(root / phase / 'uxm', 'uxm', selected, cells)
    for baseline in baselines:
        wanted = [s for s in selected if s in pd.read_csv(Path(baseline['uxm']) / 'results/lineage_fractions.csv', index_col=0).columns]
        if not wanted:
            continue
        saved = checked(pd.read_csv(Path(baseline['uxm']) / 'results/lineage_fractions.csv', index_col=0).T.loc[wanted], wanted, cells)
        delta = (native.loc[wanted] - saved) * 100
        differences.extend(dict(method='uxm', cohort=baseline['cohort'], threshold='', sample=s, cell=c, difference_pp=delta.loc[s, c]) for s in wanted for c in cells)
    for job in read_json(root / 'inputs' / f'{phase}.jobs.json'):
        native = native_fractions(root / phase / 'epidish' / job['id'], 'epidish', job['samples'], cells)
        saved = pd.read_csv(Path(job['folder']) / 'epidish_rpc_fractions.tsv', sep='\t', index_col=0).T.loc[job['samples']]
        saved = checked(saved, job['samples'], cells)
        delta = (native - saved) * 100
        differences.extend(dict(method='epidish', cohort=job['cohort'], threshold=job['threshold'], sample=s, cell=c, difference_pp=delta.loc[s, c]) for s in job['samples'] for c in cells)
    table(root / phase / 'differences.tsv', differences, ['method', 'cohort', 'threshold', 'sample', 'cell', 'difference_pp'])
    frame = pd.DataFrame(differences)
    summary = {m: {'estimates': len(g), 'maximum_absolute_difference_pp': float(g.difference_pp.abs().max()),
                   'mean_absolute_difference_pp': float(g.difference_pp.abs().mean())} for m, g in frame.groupby('method')}
    write_json(root / phase / 'comparison.json', summary)
    return summary
