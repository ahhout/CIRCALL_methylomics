#!/usr/bin/env python3
"""Check weighted PAT overlap counting, missing calls and gzip integrity."""
import gzip
import tempfile
from pathlib import Path
import unittest

from compare_healthy_wbc import scan_pat, sha, summarize


class CandidateCountsTests(unittest.TestCase):
    def test_patterns_starting_before_targets_and_multiplicity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'sample.pat.gz'
            with gzip.open(path,'wb') as handle:
                handle.write(b'chr1\t10\tCTC.T\t2\nchr1\t12\tTCC\t3\nchr1\t11\t.CC.\t5\nchr2\t100\tCCC\t9\n')
            sample = dict(sample='fixture',pat=str(path),pat_bytes=path.stat().st_size,pat_sha256=sha(path))
            targets = [dict(chr='chr1',cpg_index_1based=i) for i in [11,12,14]]
            counts = scan_pat(sample,targets)
            self.assertEqual(counts['methylated'],[0,7,3])
            self.assertEqual(counts['coverage'],[2,10,5])
            rows = [dict(region='fixture',sample='fixture',group='Healthy WBC',
                         patient_or_donor='donor',assay='WGBS',methylated=m,coverage=c)
                    for m,c in zip(counts['methylated'],counts['coverage'])]
            summary = summarize(rows)[0]
            self.assertAlmostEqual(summary['mean_fraction'],(0+.7+.6)/3)
            self.assertFalse(summary['all_cpgs_at_least_3x'])
            rows[0]['coverage']=0
            self.assertEqual(summarize(rows)[0]['mean_fraction'],'NA')

    def test_gzip_crc_failure_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'bad.pat.gz'
            with gzip.open(path,'wb') as handle: handle.write(b'chr1\t10\tCCC\t2\n')
            data=bytearray(path.read_bytes()); data[-8]^=1; path.write_bytes(data)
            sample=dict(sample='bad',pat=str(path),pat_bytes=len(data),pat_sha256=sha(path))
            with self.assertRaises(gzip.BadGzipFile):
                scan_pat(sample,[dict(chr='chr1',cpg_index_1based=11)])


if __name__=='__main__':
    unittest.main()
