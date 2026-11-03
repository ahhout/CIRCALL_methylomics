"""PAT input contract and atlas-directed proportions; fitting stays upstream."""
import numpy as np
import pandas as pd
from common import *
from fragments import write_matrices

HEME = ['Blood-B', 'Blood-T', 'Blood-NK', 'Blood-Mono+Macro', 'Blood-Granul',
        'Eryth-prog', 'Megakaryocytes']
COORDS = ['chr', 'start', 'end', 'startCpG', 'endCpG']


def profile_atlas(profile):
    data, cells = atlas_rows()
    atlas = pd.DataFrame(data)
    if profile == 'hematopoietic':
        cells = HEME
        atlas = atlas.loc[atlas.target.isin(cells)].copy()
    elif profile != 'full':
        raise ValueError(f'Unknown profile: {profile}')
    return atlas[list(atlas.columns[:8]) + cells], cells


def read_plan(path):
    path = Path(path).resolve(strict=True)
    plan = read_json(path)
    if set(plan) != {'samples', 'profiles', 'cohorts', 'fragment_rlen'}:
        raise ValueError('Plan needs samples, profiles, cohorts and fragment_rlen only')
    if plan['fragment_rlen'] != 4:
        raise ValueError('The fixed l4 atlas requires fragment_rlen=4')
    if (not plan['profiles'] or len(set(plan['profiles'])) != len(plan['profiles']) or
            not set(plan['profiles']) <= {'full', 'hematopoietic'}):
        raise ValueError('Select full and/or hematopoietic explicitly, without duplicates')
    sheet = (path.parent / plan['samples']).resolve(strict=True)
    records = rows(sheet)
    if not records or not {'sample', 'cohort', 'pat', 'csi'} <= records[0].keys():
        raise ValueError('PAT samplesheet needs sample, cohort, pat and csi')
    if len({r['sample'] for r in records}) != len(records):
        raise ValueError('Duplicate sample identifiers')
    if set(plan['cohorts']) != {r['cohort'] for r in records}:
        raise ValueError('Specify thresholds for exactly the samplesheet cohorts')
    for thresholds in plan['cohorts'].values():
        if (not isinstance(thresholds, list) or not thresholds or
                len(set(thresholds)) != len(thresholds) or
                any(type(t) is not int or t < 1 for t in thresholds)):
            raise ValueError('Coverage thresholds must be explicit unique positive integers')
    artifacts = {}
    for row in records:
        if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', row[k]) for k in ('sample', 'cohort')):
            raise ValueError('Unsafe or empty sample/cohort identifier')
        pat = (sheet.parent / row['pat']).resolve(strict=True)
        csi = (sheet.parent / row['csi']).resolve(strict=True)
        if not str(pat).endswith('.pat.gz') or csi != Path(str(pat) + '.csi').resolve(strict=True):
            raise ValueError(f'PAT/index pair does not match: {row["sample"]}')
        artifacts[row['sample']] = {'pat': identity(pat), 'csi': identity(csi, digest=True)}
    return dict(plan_path=str(path), plan_sha256=sha(path), samplesheet=str(sheet),
                samplesheet_sha256=sha(sheet), settings=plan, records=records,
                pats=artifacts, atlas_sha256=ATLAS_SHA)


def prepare_pat_inputs(root, inventory):
    dest = root / 'inputs'
    dest.mkdir()
    records = inventory['records']
    for sample, paths in inventory['pats'].items():
        for key, suffix in [('pat', '.pat.gz'), ('csi', '.pat.gz.csi')]:
            (dest / (sample + suffix)).symlink_to(paths[key]['path'])
    table(dest / 'cohort.pats.csv', [dict(sample=r['sample'], pat=str(dest / (r['sample'] + '.pat.gz')),
          csi=str(dest / (r['sample'] + '.pat.gz.csi'))) for r in records], ['sample', 'pat', 'csi'], ',')
    atlas, _ = profile_atlas('full')
    atlas[COORDS].to_csv(dest / 'markers.tsv', sep='\t', header=False, index=False)
    for profile in inventory['settings']['profiles']:
        subset, _ = profile_atlas(profile)
        target = dest / f'Atlas.U25.l4.hg38.{profile}.tsv'
        if profile == 'full':
            target.symlink_to(ATLAS.resolve(strict=True))
        else:
            subset.to_csv(target, sep='\t', index=False)
    table(root / 'samples.tsv', records, list(records[0]))
    write_json(root / 'plan.json', inventory['settings'])


def read_counts(path, atlas):
    frame = pd.read_csv(path, sep='\t', comment='#', header=None)
    if frame.shape != (len(atlas), 8):
        raise ValueError(f'Expected one eight-column count row per atlas marker: {path}')
    frame.columns = COORDS + ['U', 'X', 'M']
    expected = atlas[COORDS].copy()
    for col in COORDS[1:]:
        expected[col] = pd.to_numeric(expected[col], errors='raise')
    if not frame[COORDS].reset_index(drop=True).equals(expected.reset_index(drop=True)):
        raise ValueError(f'Count coordinates/order differ from the original atlas: {path}')
    counts = frame[['U', 'X', 'M']].to_numpy(dtype=float)
    if not np.isfinite(counts).all() or (counts < 0).any() or (counts != np.floor(counts)).any():
        raise ValueError(f'Invalid fragment counts: {path}')
    return counts


def regional_values(counts, directions):
    if counts.ndim != 2 or counts.shape != (len(directions), 3) or not set(directions) <= {'U', 'M'}:
        raise ValueError('Expected U/X/M counts and matching U/M directions')
    totals = counts.sum(axis=1)
    numerator = np.where(np.asarray(directions) == 'U', counts[:, 0], counts[:, 2])
    values = np.divide(numerator, totals, out=np.full(len(totals), np.nan), where=totals > 0)
    return totals, values


def build_mask(reference, coverage, scope, threshold):
    if type(threshold) is not int or threshold < 1:
        raise ValueError('Coverage threshold must be a positive integer')
    return (scope & np.isfinite(reference).all(axis=1) & ((reference >= 0) & (reference <= 1)).all(axis=1)
            & (coverage >= threshold).all(axis=1))


def prepared_hashes(root):
    return {str(p.relative_to(root)): sha(p) for p in (root / 'inputs').rglob('*')
            if p.is_file() and not p.is_symlink()}


def build_matrices(root, inventory):
    atlas, full_cells = profile_atlas('full')
    names = [r['sample'] for r in inventory['records']]
    coverage, values, count_hashes = [], [], {}
    for name in names:
        path = root / 'cohort/counts/fragment_counts' / f'{name}.uxm.bed.gz'
        counts = read_counts(path, atlas)
        total, value = regional_values(counts, atlas.direction)
        coverage.append(total)
        values.append(value)
        count_hashes[str(path.relative_to(root))] = sha(path)
    coverage, values = np.asarray(coverage).T, np.asarray(values).T
    reference = atlas[full_cells].apply(pd.to_numeric, errors='coerce')
    jobs = []
    for profile in inventory['settings']['profiles']:
        _, cells = profile_atlas(profile)
        ref = reference[cells].to_numpy()
        scope = np.ones(len(atlas), dtype=bool) if profile == 'full' else atlas.target.isin(cells).to_numpy()
        for cohort, thresholds in inventory['settings']['cohorts'].items():
            selected = [i for i, r in enumerate(inventory['records']) if r['cohort'] == cohort]
            wanted = [names[i] for i in selected]
            for threshold in thresholds:
                keep = build_mask(ref, coverage[:, selected], scope, threshold)
                if keep.sum() <= len(cells):
                    raise ValueError(f'Insufficient markers: {profile}/{cohort}/coverage-ge{threshold}')
                tag = f'{profile}_{cohort}_ge{threshold}'
                dest = root / 'inputs/prepared' / tag
                dest.mkdir(parents=True)
                mask = atlas.loc[keep, atlas.columns[:8]].copy()
                ids = mask['name'] + '|' + mask.direction
                mask['region_id'] = ids
                mask['minimum_cohort_coverage'] = coverage[keep][:, selected].min(axis=1)
                for i in selected:
                    mask['coverage_' + names[i]] = coverage[keep, i]
                sample = pd.DataFrame(values[keep][:, selected], index=ids, columns=wanted)
                selected_ref = pd.DataFrame(ref[keep], index=ids, columns=cells)
                if not np.isfinite(sample.values).all():
                    raise ValueError('Coverage mask retained missing sample values')
                mask.to_csv(dest / 'marker_mask.tsv', sep='\t', index=False)
                sample.to_csv(dest / 'sample_regional_values.tsv.gz', sep='\t', index_label='region_id')
                selected_ref.to_csv(dest / 'reference_matrix.tsv', sep='\t', index_label='region_id')
                paths = write_matrices(dest, tag, mask, selected_ref, sample)
                jobs.append(dict(id=tag, profile=profile, cohort=cohort, threshold=threshold,
                                 cells=cells, samples=wanted, markers=int(keep.sum()), folder=str(dest), **paths))
    write_json(root / 'inputs/cohort.jobs.json', jobs)
    write_json(root / 'counts.sha256.json', count_hashes)
    write_json(root / 'prepared.sha256.json', prepared_hashes(root))
    return jobs
