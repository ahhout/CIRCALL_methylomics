#!/usr/bin/env python3
"""Small fixtures for the actual input validator, independent of R."""
import copy
import csv
import gzip
import tempfile
import unittest
from pathlib import Path
from prepare import scan_report, validate_samples


class InputFixtures(unittest.TestCase):
    def test_manifest_rejections(self):
        with Path(__file__).with_name('config').joinpath('samples.tsv').open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        validate_samples(rows)
        for field, value in [('sample', rows[1]['sample']), ('patient', rows[1]['patient']+'wrong'),
                             ('compartment', 'cfDNA'), ('reference_name', 'GRCh37'),
                             ('report', rows[1]['report']), ('upstream_run', 'wrong')]:
            bad = copy.deepcopy(rows); bad[0][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_samples(bad)

    def test_dyad_interval_and_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp)/'input.gz'; dest = Path(tmp)/'selected.gz'
            with gzip.open(src, 'wt') as f:
                f.write('chr22\t100\t+\t0\t3\tCG\tCGA\nchr22\t101\t-\t2\t1\tCG\tCGT\n')
            totals, selected = scan_report(src, ['chr22',100,100], dest)
            self.assertEqual(selected[('chr22',100)], [2,6])
            self.assertEqual(totals['chr22'],[2,2,6])
            with gzip.open(dest,'rt') as f:
                self.assertEqual(len(f.readlines()),2)
            with gzip.open(src,'wt') as f:
                f.write('chr22\t100\t+\t-1\t3\tCG\tCGA\n')
            with self.assertRaises(ValueError):
                scan_report(src,None,dest)


if __name__ == '__main__':
    unittest.main()
