"""Corrected diagnostic restart selection and provenance failure boundaries."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import celfie_score_diagnostics as scoring
from tools.export_celfie_observed_selection import export


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def fixture(directory, best_converged=True, tie=False):
    root = Path(directory) / '2026-10-27_12-00-celfie-selection-test-v1'
    root.mkdir()
    native, counts = root / 'native.py', root / 'counts.tsv'
    native.write_text('Synthetic fixture source identity; never executed.\n')
    counts.write_text('Synthetic count source identity.\n')
    arrays = dict(alpha=np.array([[[.8, .2]], [[.25, .75]], [[.8, .2] if tie else [.75, .25]]]),
                  gamma=np.array([[[.8], [.2]]] * 3),
                  x=np.array([[7.]]), x_depths=np.array([[10.]]),
                  y=np.array([[8.], [2.]]), y_depths=np.array([[10.], [10.]]),
                  seeds=np.array([11, 22, 33]), native_scores=np.array([1., 100., 2.]))
    restarts, scores = [], []
    for index in range(3):
        score = scoring.observed_log_likelihood(
            arrays['x'], arrays['x_depths'], arrays['y'], arrays['y_depths'],
            arrays['alpha'][index], arrays['gamma'][index])
        passed = best_converged if index == 0 else True
        difference = .00001 if passed else .001
        restarts.append(dict(restart=index, seed=int(arrays['seeds'][index]), iterations=100,
                             alpha=arrays['alpha'][index].tolist(),
                             native_score=float(arrays['native_scores'][index]),
                             observed_log_likelihood=score['total'],
                             sample_log_likelihood=score['sample'], reference_log_likelihood=score['reference'],
                             converged=passed, convergence_threshold=.0001,
                             alpha_relative_change=difference, gamma_relative_change=0.,
                             final_relative_difference=difference))
        scores.append(score['total'])
    arrays['observed_scores'] = np.array(scores)
    states = root / 'scores.states.npz'
    np.savez_compressed(states, **arrays)
    record = dict(status='complete', samples=['sample_one'], cells=['Cell_A', 'Cell-B'],
                  effective_reference_counts_identical=True, effective_sample_counts_identical=True,
                  native_script=scoring.identity(native), input=scoring.identity(counts),
                  observer_script=scoring.identity(scoring.__file__), states=scoring.identity(states),
                  restarts=restarts, best_observed_restart=0, best_native_restart=1,
                  best_converged_observed_restart=0 if best_converged else 2,
                  best_observed_start_converged=best_converged, native_best_start_converged=True,
                  all_starts_converged=best_converged, observed_best_is_nonconverged=not best_converged)
    diagnostics = root / 'scores.json'
    write_json(diagnostics, record)
    return root, diagnostics, record, arrays


class ObservedSelection(unittest.TestCase):
    def setUp(self):
        guard = patch('tools.export_celfie_observed_selection.output_guard.require_clean_source', return_value='a' * 40)
        guard.start(); self.addCleanup(guard.stop)

    def test_exports_observed_winner_with_labels_and_diagnostic_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root, diagnostics, record, arrays = fixture(directory)
            preserved = root / '1_tissue_proportions.txt'
            preserved.write_text('Native winner is preserved.\n')
            before = diagnostics.read_bytes()
            result = export(diagnostics, root / '2026-11-02_13-00-selection-v1')
            self.assertEqual(result['selected_restart'], 0)
            self.assertEqual(result['native_selected_restart'], 1)
            self.assertEqual(result['selected_seed'], 11)
            self.assertTrue(result['diagnostic_only'])
            self.assertFalse(result['acceptance_granted'])
            self.assertFalse(result['biological_accuracy_validated'])
            self.assertEqual(result['reference_identity_status'], 'unresolved')
            self.assertEqual(result['status'], 'complete')
            self.assertIn('started_utc', result)
            self.assertIn('completed_utc', result)
            self.assertEqual(preserved.read_text(), 'Native winner is preserved.\n')
            self.assertEqual(diagnostics.read_bytes(), before)
            for name, expected in (
                    ('tissue_proportions.tsv', [['sample', 'Cell_A', 'Cell-B'], ['sample_one', '0.8', '0.2']]),
                    ('methylation_proportions.tsv', [['region_index', 'Cell_A', 'Cell-B'], ['0', '0.8', '0.2']])):
                with (root / '2026-11-02_13-00-selection-v1' / name).open() as handle:
                    self.assertEqual(list(csv.reader(handle, delimiter='\t')), expected)
            self.assertEqual(sorted(path.name for path in (root / '2026-11-02_13-00-selection-v1').iterdir()),
                             ['.run-ownership.json', 'export.json', 'methylation_proportions.tsv', 'tissue_proportions.tsv'])

    def test_nonconverged_best_rejected_despite_converged_alternative(self):
        with tempfile.TemporaryDirectory() as directory:
            root, diagnostics, _, _ = fixture(directory, best_converged=False)
            with self.assertRaisesRegex(ValueError, 'no converged fallback'):
                export(diagnostics, root / '2026-11-02_13-00-selection-v1')
            result = json.loads((root / '2026-11-02_13-00-selection-v1/export.json').read_text())
            self.assertEqual(result['status'], 'failed')
            self.assertIn('ended_utc', result)
            self.assertFalse((root / '2026-11-02_13-00-selection-v1/tissue_proportions.tsv').exists())

    def test_tie_selects_first_restart_deterministically(self):
        with tempfile.TemporaryDirectory() as directory:
            root, diagnostics, _, _ = fixture(directory, tie=True)
            first = export(diagnostics, root / '2026-11-02_13-00-selection-one-v1')
            second = export(diagnostics, root / '2026-11-02_13-00-selection-two-v1')
            self.assertEqual(first['selected_restart'], 0)
            self.assertEqual(second['selected_restart'], 0)
            self.assertEqual(first['outputs']['fractions']['sha256'], second['outputs']['fractions']['sha256'])

    def test_tampered_npz_rejected_by_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            root, diagnostics, record, _ = fixture(directory)
            states = Path(record['states']['path'])
            states.write_bytes(states.read_bytes() + b'tampered')
            with self.assertRaisesRegex(ValueError, 'states identity'):
                export(diagnostics, root / '2026-11-02_13-00-selection-v1')

    def test_changed_arrays_rejected_even_with_updated_archive_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            root, diagnostics, record, arrays = fixture(directory)
            arrays['gamma'][0, 0, 0] = .6
            np.savez_compressed(record['states']['path'], **arrays)
            record['states'] = scoring.identity(record['states']['path'])
            write_json(diagnostics, record)
            with self.assertRaisesRegex(ValueError, 'Stored observed score'):
                export(diagnostics, root / '2026-11-02_13-00-selection-v1')

    def test_tampered_ranking_or_source_rejected(self):
        for what in ('ranking', 'source'):
            with self.subTest(what=what), tempfile.TemporaryDirectory() as directory:
                root, diagnostics, record, _ = fixture(directory)
                if what == 'ranking':
                    record['best_observed_restart'] = 1
                    write_json(diagnostics, record)
                else:
                    Path(record['native_script']['path']).write_text('changed source')
                with self.assertRaisesRegex(ValueError, 'ranking|native_script identity'):
                    export(diagnostics, root / '2026-11-02_13-00-selection-v1')

    def test_nonfinite_score_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root, diagnostics, record, arrays = fixture(directory)
            arrays['native_scores'][1] = np.nan
            np.savez_compressed(record['states']['path'], **arrays)
            record['states'] = scoring.identity(record['states']['path'])
            write_json(diagnostics, record)
            with self.assertRaisesRegex(ValueError, 'scores must all be finite'):
                export(diagnostics, root / '2026-11-02_13-00-selection-v1')

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            root, diagnostics, _, _ = fixture(directory)
            output = root / '2026-11-02_13-00-selection-v1'
            output.mkdir()
            (output / 'existing.txt').write_text('preserved')
            with self.assertRaises(FileExistsError):
                export(diagnostics, output)
            self.assertEqual((output / 'existing.txt').read_text(), 'preserved')


if __name__ == '__main__':
    unittest.main()
