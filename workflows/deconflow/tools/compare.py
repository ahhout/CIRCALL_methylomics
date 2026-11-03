#!/usr/bin/env python3
"""Compare saved native DecoNFlow estimates with immutable historical fits; never refit."""
import argparse
import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import *
from native import checked, native_fractions



def differences(native, baseline, method, label):
    records = []
    for sample in baseline.index:
        for cell in baseline.columns:
            a, b = float(native.loc[sample, cell]), float(baseline.loc[sample, cell])
            records.append(dict(method=method, baseline=label, sample=sample, cell_type=cell,
                                native_fraction=a, baseline_fraction=b,
                                difference_pp=100 * (a - b), absolute_difference_pp=100 * abs(a - b)))
    return records


def compare(root, analysis, phase):
    if not (analysis / (phase.upper() + '_COMPLETE')).is_file():
        raise ValueError('Native phase has not passed both methods')
    from analysis import selected_samples
    metadata = selected_samples(rows(analysis / 'samples.tsv'), phase)
    samples = [r['sample'] for r in metadata]
    _, cells = atlas_rows()
    native = {method: native_fractions(analysis / phase / method, method, samples, cells) for method in ('uxm', 'epidish')}
    records, provenance, support, warnings = [], [], [], []
    seen = set()
    for cohort in read_json(analysis / 'baselines.json'):
        uxm = Path(cohort['uxm']) / 'results/lineage_fractions.csv'
        old = pd.read_csv(uxm, index_col=0).T
        selected = [s for s in samples if s in old.index]
        if not selected:
            continue
        if seen.intersection(selected):
            raise ValueError('Overlapping baseline cohorts')
        seen.update(selected)
        old = checked(old, list(old.index), cells).loc[selected]
        records += differences(native['uxm'], old, 'uxm', cohort['cohort'] + ':native')
        input_manifest = Path(cohort['epidish']) / 'provenance/input_manifest.json'
        saved = read_json(input_manifest)
        if sha(uxm) != saved['full_uxm_run']['native_result_sha256']:
            raise ValueError('Saved native UXM checksum disagrees with baseline provenance')
        provenance += [identity(uxm, True), identity(input_manifest, True)]
        folders = sorted((Path(cohort['epidish']) / 'comparisons/full').glob('coverage-ge*'))
        if not folders:
            raise ValueError('No saved full-atlas EpiDISH thresholds')
        for folder in folders:
            manifest = read_json(folder / 'manifest.json')
            if manifest['profile'] != 'full':
                raise ValueError('Unexpected restricted baseline profile')
            for filename in ('epidish_rpc_fractions.tsv', 'epidish_warnings.tsv', 'marker_mask.tsv'):
                expected = next(r['sha256'] for r in manifest['files'] if r['name'] == filename)
                if sha(folder / filename) != expected:
                    raise ValueError(f'Saved baseline artifact changed: {folder / filename}')
                provenance.append(identity(folder / filename, True))
            old = pd.read_csv(folder / 'epidish_rpc_fractions.tsv', sep='\t', index_col=0).T
            old = checked(old, manifest['samples'], cells).loc[selected]
            label = f'{cohort["cohort"]}:coverage-ge{manifest["coverage_threshold"]}'
            records += differences(native['epidish'], old, 'epidish', label)
            for sample in selected:
                support.append({'source': label, 'method': 'epidish_baseline', 'sample': sample,
                                'markers': manifest['marker_count'], 'rpc_maxit': saved.get('epidish_maxit', 'not-recorded')})
            warnings += [dict(r, source=label) for r in rows(folder / 'epidish_warnings.tsv')]
    if seen != set(samples):
        raise ValueError('Some native samples have no saved baseline')
    df = pd.DataFrame(records)
    df.to_csv(root / 'all_differences.tsv', sep='\t', index=False, float_format='%.12g')
    summary = df.groupby(['method', 'baseline', 'sample'], sort=False).agg(
        mean_absolute_difference_pp=('absolute_difference_pp', 'mean'),
        max_absolute_difference_pp=('absolute_difference_pp', 'max'),
        sum_absolute_difference_pp=('absolute_difference_pp', 'sum'))
    summary.to_csv(root / 'difference_summary.tsv', sep='\t', float_format='%.12g')
    # Native coverage support is measured from the actual matrix, without changing masks.
    mixes = list((analysis / phase / 'epidish/preprocessing/bulk_samples_preprocessing').glob('*sample_mix.csv'))
    for path in mixes:
        frame = pd.read_csv(path, index_col=0)
        for sample in samples:
            if sample in frame.columns:
                support.append({'source': str(path), 'method': 'epidish_native', 'sample': sample,
                                'markers': int(frame[sample].notna().sum()), 'rpc_maxit': 10000})
    table(root / 'marker_support.tsv', support, ['source', 'method', 'sample', 'markers', 'rpc_maxit'])
    write_json(root / 'baseline_warnings.json', warnings)
    # Keep native warning lines with exact originating log/task paths.
    native_logs, native_uxm_counts = [], []
    for method in ('uxm', 'epidish'):
        trace = analysis / phase / method / 'pipeline_info' / f'execution_trace_{phase}.txt'
        provenance.append(identity(trace, True))
        for row in rows(trace):
            workdir = row.get('workdir')
            if workdir:
                if method == 'uxm' and 'UXM_DECONV' in row['name']:
                    for homog in sorted((Path(workdir) / '.tmp').glob('l*/*.homog.gz')):
                        sample = homog.name.split('.mem.')[0]
                        if sample not in samples:
                            raise ValueError('Unrecognized sample in native UXM support table')
                        counts = pd.read_csv(homog, sep='\t', header=None, comment='#',
                                             names=['chr', 'start', 'end', 'startCpG', 'endCpG', 'U', 'X', 'M'])
                        counts.insert(0, 'sample', sample)
                        native_uxm_counts.append(counts)
                        support.append({'source': str(homog), 'method': 'uxm_native', 'sample': sample,
                                        'markers': int((counts[['U', 'X', 'M']].sum(axis=1) > 0).sum()),
                                        'rpc_maxit': 'not-applicable'})
                for filename in ('.command.err', '.command.out'):
                    p = Path(workdir) / filename
                    if p.is_file():
                        native_logs += [{'file': str(p), 'line': line} for line in p.read_text(errors='replace').splitlines()
                                        if re.search('warning|converg|error', line, re.I)]
    write_json(root / 'native_warnings.json', native_logs)
    table(root / 'marker_support.tsv', support, ['source', 'method', 'sample', 'markers', 'rpc_maxit'])
    if native_uxm_counts:
        pd.concat(native_uxm_counts).to_csv(root / 'native_uxm_marker_counts.tsv', sep='\t', index=False)
    os.environ.setdefault('MPLCONFIGDIR', str(root / '.plot-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(root / 'scatterplots.pdf') as pdf:
        for (method, baseline), frame in df.groupby(['method', 'baseline'], sort=False):
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(100 * frame.baseline_fraction, 100 * frame.native_fraction, s=10, alpha=.5)
            ax.plot([0, 100], [0, 100], color='grey', lw=.7)
            ax.set(xlabel='Saved fraction (%)', ylabel='Native DecoNFlow fraction (%)', title=f'{method}: {baseline}')
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    with PdfPages(root / 'compositions.pdf') as pdf:
        for method, frame in native.items():
            fig, ax = plt.subplots(figsize=(15, 8))
            frame.plot.bar(stacked=True, ax=ax, colormap='gist_ncar')
            ax.set(ylabel='Fraction', title=f'Native DecoNFlow {method}: descriptive normal-cell composition')
            ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=6)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    write_json(root / 'comparison.provenance.json', {'analysis': str(analysis), 'phase': phase, 'sources': provenance,
               'native_fingerprint': identity(analysis / 'input_fingerprint.json', True),
               'difference_definition': '100 * (native - saved)', 'discrepancy_cutoff': None,
               'interpretation': 'EpiDISH comparisons combine representation, atlas construction, marker handling and software differences. PB is diagnostic material; these are not tumour-specific fractions or DMR covariates.'})
    finish(root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--phase', choices=['triplet', 'cohort'], default='cohort')
    args = parser.parse_args()
    root = new_run(args.output)
    try:
        compare(root, args.analysis.resolve(strict=True), args.phase)
    except Exception as exc:
        finish(root, 'failed', error=str(exc))
        raise


if __name__ == '__main__':
    main()
