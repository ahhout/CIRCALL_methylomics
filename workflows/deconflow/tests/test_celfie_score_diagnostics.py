"""Independent likelihood checks and unchanged-native observer boundaries."""
import json
import math
from pathlib import Path
import runpy
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.celfie_score_diagnostics import observed_log_likelihood, observe


def scalar_oracle(x, depths, y, reference_depths, alpha, gamma):
    """Explicit scalar Bernoulli sums, independent of the vectorized scorer."""
    sample, reference = 0.0, 0.0
    for individual in range(len(x)):
        for site in range(len(x[individual])):
            probability = sum(alpha[individual][cell] * gamma[cell][site]
                              for cell in range(len(y)))
            if x[individual][site]:
                sample += x[individual][site] * math.log(probability)
            if depths[individual][site] - x[individual][site]:
                sample += (depths[individual][site] - x[individual][site]) * math.log(1 - probability)
    for cell in range(len(y)):
        for site in range(len(y[cell])):
            if y[cell][site]:
                reference += y[cell][site] * math.log(gamma[cell][site])
            if reference_depths[cell][site] - y[cell][site]:
                reference += (reference_depths[cell][site] - y[cell][site]) * math.log(1 - gamma[cell][site])
    return dict(total=sample + reference, sample=sample, reference=reference)


class Likelihood(unittest.TestCase):
    def setUp(self):
        self.x = [[7, 2], [3, 5]]
        self.depths = [[10, 8], [6, 9]]
        self.y = [[8, 1], [2, 7]]
        self.reference_depths = [[10, 10], [10, 10]]
        self.alpha = [[.8, .2], [.3, .7]]
        self.gamma = [[.7, .2], [.25, .75]]

    def test_independent_scalar_likelihood_oracle(self):
        arguments = (self.x, self.depths, self.y, self.reference_depths, self.alpha, self.gamma)
        expected = scalar_oracle(*arguments)
        actual = observed_log_likelihood(*arguments)
        for name in expected:
            self.assertAlmostEqual(actual[name], expected[name], places=12)

    def test_sample_replication_does_not_repeat_reference(self):
        original = observed_log_likelihood(self.x, self.depths, self.y,
                                          self.reference_depths, self.alpha, self.gamma)
        repeated = observed_log_likelihood(self.x * 3, self.depths * 3, self.y,
                                          self.reference_depths, self.alpha * 3, self.gamma)
        self.assertAlmostEqual(repeated['sample'], 3 * original['sample'])
        self.assertEqual(repeated['reference'], original['reference'])
        self.assertAlmostEqual(repeated['total'], 3 * original['sample'] + original['reference'])

    def test_zero_times_log_zero_and_impossible_observations(self):
        with np.errstate(all='raise'):
            result = observed_log_likelihood([[0, 2]], [[0, 2]], [[0, 2]],
                                             [[0, 2]], [[1]], [[0, 1]])
            self.assertEqual(result, dict(total=0, sample=0, reference=0))
            impossible = observed_log_likelihood([[1]], [[1]], [[0]],
                                                 [[0]], [[1]], [[0]])
            self.assertEqual(impossible['sample'], -math.inf)
            self.assertEqual(impossible['total'], -math.inf)
            self.assertEqual(impossible['reference'], 0)
            impossible_reference = observed_log_likelihood([[0]], [[0]], [[1]],
                                                           [[1]], [[1]], [[0]])
            self.assertEqual(impossible_reference['reference'], -math.inf)

    def test_equal_profiles_make_observed_likelihood_alpha_invariant(self):
        arguments = ([[5]], [[10]], [[5], [5]], [[10], [10]])
        equal = observed_log_likelihood(*arguments, [[.5, .5]], [[.5], [.5]])
        concentrated = observed_log_likelihood(*arguments, [[.99, .01]], [[.5], [.5]])
        self.assertEqual(equal, concentrated)
        self.assertAlmostEqual(equal['total'], 30 * math.log(.5))

    def test_invalid_counts_or_parameters_rejected(self):
        with self.assertRaisesRegex(ValueError, 'sample'):
            observed_log_likelihood([[2]], [[1]], [[1]], [[2]], [[1]], [[.5]])
        with self.assertRaisesRegex(ValueError, 'alpha/gamma'):
            observed_log_likelihood([[1]], [[2]], [[1]], [[2]], [[.8]], [[.5]])


NATIVE = '''import json
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
args = SimpleNamespace(random_restarts=3, input_path=sys.argv[1])
samples, tissues = ['sample'], ['A', 'B']
initial_seed = 42
x, x_depths = np.array([[7.]]), np.array([[10.]])
y, y_depths = np.array([[8.], [2.]]), np.array([[10.], [10.]])
def em(x, x_depths, y, y_depths, seed):
    i = 9
    convergence_criteria = .0001
    fractions, score, difference = {11: ([.8, .2], 1., .001),
                                   22: ([.25, .75], 100., .00001),
                                   33: ([.75, .25], 2., .00002)}[seed]
    alpha, gamma = np.array([fractions]), y / y_depths
    alpha_diff, gamma_diff = difference, 0.
    a = alpha + np.array([[.001, -.001]]) if difference < convergence_criteria else alpha.copy()
    g = gamma.copy()
    return alpha, gamma, score
random_restarts = []
for seed in (11, 22, 33):
    alpha, gamma, ll = em(x, x_depths, y, y_depths, seed)
    random_restarts.append((ll, alpha, gamma))
ll_max, alpha_max, gamma_max = max(random_restarts, key=lambda result: result[0])
output_alpha_file = str(Path(sys.argv[2]) / 'alpha.txt')
output_gamma_file = str(Path(sys.argv[2]) / 'gamma.txt')
Path(output_alpha_file).write_text(json.dumps(alpha_max.tolist()))
Path(output_gamma_file).write_text(json.dumps(gamma_max.tolist()))
'''


class Observer(unittest.TestCase):
    def prepare(self, directory, source=NATIVE):
        root = Path(directory) / '2026-10-27_12-00-celfie-score-test-v1'
        root.mkdir()
        script, counts = root / 'native.py', root / 'counts.tsv'
        script.write_text(source)
        counts.write_text('synthetic count input\n')
        return root, script, counts

    def test_preserves_native_outputs_and_separates_diagnostic_rankings(self):
        with tempfile.TemporaryDirectory() as directory:
            root, script, counts = self.prepare(directory)
            plain, observed = root / 'plain', root / 'observed'
            plain.mkdir()
            observed.mkdir()
            original_argv, original_profile, errors = sys.argv, sys.getprofile(), np.geterr()
            try:
                sys.argv = [str(script), str(counts), str(plain), '1']
                runpy.run_path(str(script), run_name='__main__')
            finally:
                sys.argv = original_argv
                np.seterr(**errors)
            source_before = script.read_bytes()
            result = observe(script, root / 'scores.json', [counts, observed, '1'])
            self.assertIs(sys.argv, original_argv)
            self.assertIs(sys.getprofile(), original_profile)
            self.assertEqual(np.geterr(), errors)
            self.assertEqual(script.read_bytes(), source_before)
            for name in ('alpha.txt', 'gamma.txt'):
                self.assertEqual((plain / name).read_bytes(), (observed / name).read_bytes())
            self.assertEqual(result['best_native_restart'], 1)
            self.assertEqual(result['best_observed_restart'], 0)
            self.assertEqual(result['best_converged_observed_restart'], 2)
            self.assertFalse(result['best_observed_start_converged'])
            self.assertTrue(result['observed_best_is_nonconverged'])
            self.assertFalse(result['all_starts_converged'])
            self.assertFalse(result['acceptance_granted'])
            self.assertTrue(result['effective_reference_counts_identical'])
            self.assertEqual(result['status'], 'complete')
            self.assertIn('started_utc', result)
            self.assertIn('completed_utc', result)
            self.assertEqual(result['restarts'][0]['candidate_alpha_max_absolute_difference_from_returned'], 0)
            self.assertAlmostEqual(result['restarts'][1]['candidate_alpha_max_absolute_difference_from_returned'], .001)
            self.assertEqual(len(result['input']['sha256']), 64)
            self.assertEqual(json.loads((root / 'scores.json').read_text()), result)
            with np.load(result['states']['path'], allow_pickle=False) as states:
                self.assertEqual(states['alpha'].shape, (3, 1, 2))
                self.assertEqual(states['gamma'].shape, (3, 2, 1))
                np.testing.assert_array_equal(states['seeds'], [11, 22, 33])
                np.testing.assert_array_equal(states['y'], [[8.], [2.]])
            self.assertFalse(any(root.glob('*VALIDATED*')))

    def test_no_converged_start_has_no_converged_best(self):
        source = NATIVE.replace('convergence_criteria = .0001', 'convergence_criteria = 1e-10')
        with tempfile.TemporaryDirectory() as directory:
            root, script, counts = self.prepare(directory, source)
            result = observe(script, root / 'scores.json', [counts, root, '1'])
            self.assertIsNone(result['best_converged_observed_restart'])
            self.assertTrue(result['observed_best_is_nonconverged'])

    def test_different_effective_reference_counts_prevent_score_comparison(self):
        source = NATIVE.replace('    i = 9', '    y[0, 0] += .1\n    i = 9')
        with tempfile.TemporaryDirectory() as directory:
            root, script, counts = self.prepare(directory, source)
            with self.assertRaisesRegex(ValueError, 'reference counts differ'):
                observe(script, root / 'scores.json', [counts, root, '1'])
            record = json.loads((root / 'scores.json').read_text())
            self.assertFalse(record['effective_reference_counts_identical'])
            self.assertEqual(record['status'], 'failed')
            self.assertNotIn('best_observed_restart', record)
            self.assertFalse((root / 'scores.states.npz').exists())
            self.assertTrue((root / 'alpha.txt').exists())

    def test_native_exception_survives_and_profile_is_restored(self):
        with tempfile.TemporaryDirectory() as directory:
            root, script, counts = self.prepare(directory, "def em():\n    raise RuntimeError('native failure')\nem()\n")
            before = sys.getprofile()
            with self.assertRaisesRegex(RuntimeError, 'native failure'):
                observe(script, root / 'scores.json', [counts, root, '1'])
            self.assertIs(sys.getprofile(), before)
            record = json.loads((root / 'scores.json').read_text())
            self.assertEqual(record['status'], 'failed')
            self.assertIn('started_utc', record)
            self.assertIn('ended_utc', record)

    def test_existing_outputs_and_invalid_run_date_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root, script, counts = self.prepare(directory)
            output = root / 'scores.json'
            output.write_text('existing')
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                observe(script, output, [counts, root, '1'])
            self.assertEqual(output.read_text(), 'existing')
            invalid = Path(directory) / '2026-02-30_12-00-celfie-score-test-v1'
            invalid.mkdir()
            with self.assertRaises(ValueError):
                observe(script, invalid / 'scores.json', [counts, root, '1'])

    def test_existing_native_output_is_protected_before_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            root, script, counts = self.prepare(directory)
            native_output = root / '8_tissue_proportions.txt'
            native_output.write_text('preserved native result')
            with patch.object(Path, 'cwd', return_value=root):
                with self.assertRaisesRegex(ValueError, 'overwrite native output'):
                    observe(script, root / 'scores.json', [counts, root, '1', '-p', '8'])
            self.assertEqual(native_output.read_text(), 'preserved native result')
            self.assertFalse((root / 'scores.json').exists())


if __name__ == '__main__':
    unittest.main()
