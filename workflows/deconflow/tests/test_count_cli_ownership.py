"""Direct counter CLI reserves a fresh bundle; library work uses parent ownership."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import output_guard
from tools import pat_counts


class CounterCLI(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(dir='/tmp'); self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.pat = self.root / 'fixture.pat.gz'; self.pat.write_text('fixture payload')
        self.atlas = self.root / 'atlas.tsv'; self.atlas.write_text('fixture metadata')
        self.output = self.root / '2026-11-02_13-00-count-fixture-v1'
        argv = ['pat_counts.py', '--pat', str(self.pat), '--atlas', str(self.atlas), '--output', str(self.output)]
        for patcher in (patch.object(sys, 'argv', argv),
                        patch.object(output_guard, 'require_clean_source', return_value='a' * 40)):
            patcher.start(); self.addCleanup(patcher.stop)

    def test_complete_bundle_cannot_be_reused(self):
        with patch.object(pat_counts, 'count_pat', return_value={'status': 'completed'}) as count:
            pat_counts.main()
            self.assertEqual(output_guard.read_record(self.output)['status'], 'completed')
            with self.assertRaises(FileExistsError):
                pat_counts.main()
            self.assertEqual(count.call_count, 1)

    def test_failed_initialization_retains_owned_failure(self):
        with patch.object(pat_counts, 'count_pat', side_effect=ValueError('mock failure')):
            with self.assertRaisesRegex(ValueError, 'mock failure'):
                pat_counts.main()
        self.assertEqual(output_guard.read_record(self.output)['status'], 'failed')
        with self.assertRaises(FileExistsError):
            pat_counts.main()


if __name__ == '__main__':
    unittest.main()
