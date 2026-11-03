"""Input scientific contracts that matter when there is no saved UXM run."""
import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import table, write_json
from pat_inputs import read_plan, read_counts, regional_values, build_mask, profile_atlas


class PatInputs(unittest.TestCase):
    def test_new_samples_need_no_historical_cohort_or_patient_triplets(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            records = []
            for sample in ['new_a', 'new_b']:
                pat = root / (sample + '.pat.gz')
                pat.write_bytes(b'fixture')
                Path(str(pat) + '.csi').write_bytes(b'index')
                records.append(dict(sample=sample, cohort='pilot', pat=pat.name, csi=pat.name + '.csi'))
            table(root / 'samples.tsv', records, ['sample', 'cohort', 'pat', 'csi'])
            plan = dict(samples='samples.tsv', profiles=['hematopoietic'], cohorts={'pilot': [1, 3]}, fragment_rlen=4)
            write_json(root / 'plan.json', plan)
            actual = read_plan(root / 'plan.json')
            self.assertEqual(len(actual['records']), 2)
            self.assertEqual(actual['pats']['new_a']['pat']['path'], str(root / 'new_a.pat.gz'))
            plan['cohorts'] = {}
            write_json(root / 'plan.json', plan)
            with self.assertRaises(ValueError):
                read_plan(root / 'plan.json')
            plan['cohorts'] = {'pilot': [0]}
            write_json(root / 'plan.json', plan)
            with self.assertRaises(ValueError):
                read_plan(root / 'plan.json')
            plan['cohorts'] = {'pilot': [1]}
            plan['fragment_rlen'] = 3
            write_json(root / 'plan.json', plan)
            with self.assertRaises(ValueError):
                read_plan(root / 'plan.json')

    def test_direction_and_zero_coverage(self):
        totals, values = regional_values(np.array([[2., 3., 5.], [2., 3., 5.], [0., 0., 0.]]), ['U', 'M', 'U'])
        np.testing.assert_array_equal(totals, [10, 10, 0])
        np.testing.assert_allclose(values[:2], [.2, .5])
        self.assertTrue(np.isnan(values[2]))

    def test_mask_uses_all_cohort_samples_and_selected_cell_reference(self):
        reference = np.array([[.1, .9], [.2, np.nan], [.3, .7], [.4, .6]])
        coverage = np.array([[3, 4], [4, 4], [4, 0], [3, 4]])
        scope = np.array([True, True, True, False])
        np.testing.assert_array_equal(build_mask(reference, coverage, scope, 3), [True, False, False, False])
        np.testing.assert_array_equal(build_mask(reference[:, :1], coverage, scope, 3), [True, True, False, False])

    def test_count_coordinates_and_integer_counts_are_checked(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            atlas = pd.DataFrame([['chr1', 10, 20, 2, 5]], columns=['chr', 'start', 'end', 'startCpG', 'endCpG'])
            path = root / 'sample.uxm.bed.gz'
            def dump(row):
                pd.DataFrame([row]).to_csv(path, sep='\t', index=False, header=False)
            dump(['chr1', 10, 20, 2, 5, 2, 1, 4])
            np.testing.assert_array_equal(read_counts(path, atlas), [[2, 1, 4]])
            dump(['chr1', 11, 20, 2, 5, 2, 1, 4])
            with self.assertRaises(ValueError):
                read_counts(path, atlas)
            dump(['chr1', 10, 20, 2, 5, 2.5, 1, 4])
            with self.assertRaises(ValueError):
                read_counts(path, atlas)

    def test_hematopoietic_profile_preserves_original_targets_and_seven_cells(self):
        full, _ = profile_atlas('full')
        heme, cells = profile_atlas('hematopoietic')
        self.assertEqual(len(cells), 7)
        self.assertEqual(len(heme), 175)
        pd.testing.assert_frame_equal(heme, full.loc[full.target.isin(cells), list(full.columns[:8]) + cells])


if __name__ == '__main__':
    unittest.main()
