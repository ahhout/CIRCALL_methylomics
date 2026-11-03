#!/usr/bin/env python3
"""Export a converged observed-score winner as an explicitly diagnostic result."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import celfie_score_diagnostics as scoring
from common import output_guard


def verify_identity(record, label):
    actual = scoring.identity(record['path'])
    if actual != record:
        raise ValueError(f'{label} identity differs from recorded provenance')
    return actual


def labels(values, label):
    if (not isinstance(values, list) or not values
            or any(not isinstance(value, str) or not value
                   or any(character in value for character in '\t\n\r') for value in values)
            or len(values) != len(set(values))):
        raise ValueError(f'Invalid or duplicate {label} labels')
    return values


def same_score(actual, recorded, label):
    if (recorded is None or not np.isfinite(recorded) or not np.isfinite(actual)
            or not np.isclose(actual, recorded, rtol=1e-12, atol=1e-8)):
        raise ValueError(f'{label} is nonfinite or differs from recomputed score')


def checked_states(document):
    """Verify provenance and recalculate every ranking before selecting a state."""
    if document.get('status') != 'complete':
        raise ValueError('A completed diagnostic observation is required')
    if (document.get('effective_reference_counts_identical') is not True
            or document.get('effective_sample_counts_identical') is not True):
        raise ValueError('Identical effective counts across restarts are required')
    sources = {name: verify_identity(document[name], name)
               for name in ('native_script', 'observer_script', 'input', 'states')}
    if sources['observer_script']['sha256'] != scoring.identity(scoring.__file__)['sha256']:
        raise ValueError('Installed scorer differs from the recorded observer source')
    samples, cells = labels(document['samples'], 'sample'), labels(document['cells'], 'cell')
    restarts = document['restarts']
    if not restarts:
        raise ValueError('No restart states are recorded')
    required = ('alpha', 'gamma', 'x', 'x_depths', 'y', 'y_depths',
                'seeds', 'native_scores', 'observed_scores')
    with np.load(sources['states']['path'], allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in required}
    number = len(restarts)
    if (arrays['alpha'].shape != (number, len(samples), len(cells))
            or arrays['gamma'].ndim != 3 or arrays['gamma'].shape[:2] != (number, len(cells))
            or any(arrays[name].shape != (number,)
                   for name in ('seeds', 'native_scores', 'observed_scores'))):
        raise ValueError('State arrays disagree with restart/sample/cell labels')
    if not np.isfinite(arrays['native_scores']).all():
        raise ValueError('Native restart scores must all be finite')
    scores, converged = [], []
    for index, restart in enumerate(restarts):
        if restart['restart'] != index or restart['seed'] != arrays['seeds'][index]:
            raise ValueError('Restart order or seed differs between JSON and states')
        if not np.array_equal(arrays['alpha'][index], np.asarray(restart['alpha'])):
            raise ValueError('Terminal fractions differ between JSON and states')
        actual = scoring.observed_log_likelihood(
            arrays['x'], arrays['x_depths'], arrays['y'], arrays['y_depths'],
            arrays['alpha'][index], arrays['gamma'][index])
        same_score(actual['total'], arrays['observed_scores'][index], 'Stored observed score')
        for component, field in (('total', 'observed_log_likelihood'),
                                 ('sample', 'sample_log_likelihood'),
                                 ('reference', 'reference_log_likelihood')):
            same_score(actual[component], restart[field], field)
        same_score(arrays['native_scores'][index], restart['native_score'], 'Native score')
        changes = [restart['alpha_relative_change'], restart['gamma_relative_change']]
        threshold = restart['convergence_threshold']
        if (not all(np.isfinite(value) and value >= 0 for value in changes)
                or not np.isfinite(threshold) or threshold <= 0
                or not isinstance(restart['iterations'], int) or restart['iterations'] < 1):
            raise ValueError('Invalid recorded convergence information')
        difference = sum(changes)
        same_score(difference, restart['final_relative_difference'], 'Final parameter change')
        passed = bool(difference < threshold)
        if restart['converged'] is not passed:
            raise ValueError('Convergence flag disagrees with the recorded parameter change')
        converged.append(passed)
        scores.append(actual['total'])
    # Python max returns the first index on exact ties, matching the observer.
    best = max(range(number), key=lambda index: scores[index])
    native_best = max(range(number), key=lambda index: arrays['native_scores'][index])
    passing = [index for index in range(number) if converged[index]]
    converged_best = max(passing, key=lambda index: scores[index]) if passing else None
    expected = dict(best_observed_restart=best, best_native_restart=native_best,
                    best_converged_observed_restart=converged_best,
                    best_observed_start_converged=converged[best],
                    native_best_start_converged=converged[native_best],
                    all_starts_converged=all(converged),
                    observed_best_is_nonconverged=not converged[best])
    if any(document.get(key) != value for key, value in expected.items()):
        raise ValueError('Recorded ranking or convergence summary differs from recomputed states')
    if not converged[best]:
        raise ValueError('Highest observed-score restart did not converge; no converged fallback is allowed')
    return arrays, best, sources


def export(diagnostics, output):
    revision = output_guard.require_clean_source(Path(__file__))
    started = scoring.utc()
    diagnostics = Path(diagnostics).resolve(strict=True)
    output = Path(output).absolute()
    scoring._run_root(output / 'export.json')
    if not output.parent.is_dir():
        raise ValueError('The output parent must already exist')
    output_guard.reserve_path(output, {'project_revision': revision, 'source': scoring.identity(__file__), 'diagnostics': scoring.identity(diagnostics)})
    manifest = dict(
        status='running', started_utc=started, diagnostic_only=True,
        acceptance_granted=False, biological_accuracy_validated=False,
        reference_identity_status='unresolved',
        reference_limitation='Reconstructed megakaryocyte counts remain unreconciled with the UXM atlas',
        selection='Highest finite observed-data likelihood; first restart on exact ties; winner must converge; no fallback',
        native_updates_changed=False, native_outputs_replaced=False,
        source=scoring.identity(__file__), scoring_source=scoring.identity(scoring.__file__),
        diagnostics=scoring.identity(diagnostics))
    manifest_path = output / 'export.json'

    def save():
        manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')

    save()
    try:
        document = json.loads(diagnostics.read_text())
        arrays, selected, sources = checked_states(document)
        if scoring.identity(diagnostics) != manifest['diagnostics']:
            raise ValueError('Diagnostic JSON changed during verification')
        fractions = output / 'tissue_proportions.tsv'
        methylation = output / 'methylation_proportions.tsv'
        for path, label, indices, values in (
                (fractions, 'sample', document['samples'], arrays['alpha'][selected]),
                (methylation, 'region_index', range(arrays['gamma'].shape[2]), arrays['gamma'][selected].T)):
            with path.open('x', newline='') as handle:
                writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                writer.writerow([label, *document['cells']])
                for name, row in zip(indices, values):
                    writer.writerow([name, *[repr(float(value)) for value in row]])
        manifest.update(
            status='complete', completed_utc=scoring.utc(), selected_restart=selected,
            selected_seed=int(arrays['seeds'][selected]),
            selected_observed_log_likelihood=float(arrays['observed_scores'][selected]),
            selected_start_converged=True, native_selected_restart=document['best_native_restart'],
            samples=document['samples'], cells=document['cells'], regions=arrays['gamma'].shape[2],
            gamma_row_index='Zero-based native region order; no coordinate conversion',
            verified_sources=sources,
            outputs=dict(fractions=scoring.identity(fractions), methylation=scoring.identity(methylation)))
    except BaseException as exc:
        manifest.update(status='failed', ended_utc=scoring.utc(), error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        save()
        output_guard.finish(output, manifest['status'] == 'complete')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--diagnostics', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    export(args.diagnostics, args.output)


if __name__ == '__main__':
    main()
