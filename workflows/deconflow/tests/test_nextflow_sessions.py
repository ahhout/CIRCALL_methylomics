"""Session receipts must identify the original local cache, never the latest run."""
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import save_nextflow_session, resume_nextflow_session, current_nextflow_session

class NextflowSessions(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir='/tmp')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.session = '01234567-1234-1234-1234-123456789abc'
        (self.root/'.nextflow/cache'/self.session/'db').mkdir(parents=True)
        (self.root/'.nextflow/history').write_text('fixture ' + self.session + '\n')
        (self.root/'.nextflow.log').write_text('Session UUID: ' + self.session + '\n')

    def test_exact_session_round_trip(self):
        self.assertEqual(save_nextflow_session(self.root), self.session)
        self.assertEqual(resume_nextflow_session(self.root), self.session)

    def test_missing_receipt_refuses_adoption(self):
        with self.assertRaisesRegex(ValueError, 'receipt'):
            resume_nextflow_session(self.root)

    def test_changed_or_conflicting_session_refused(self):
        save_nextflow_session(self.root)
        (self.root/'.nextflow.log').write_text('Session UUID: ' + self.session + '\nSession UUID: aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa\n')
        with self.assertRaisesRegex(ValueError, 'conflicting'):
            resume_nextflow_session(self.root)

    def test_missing_cache_refused(self):
        save_nextflow_session(self.root)
        (self.root/'.nextflow/cache'/self.session/'db').rmdir()
        with self.assertRaisesRegex(ValueError, 'cache/history'):
            resume_nextflow_session(self.root)

    def test_history_drift_refused(self):
        save_nextflow_session(self.root)
        (self.root/'.nextflow/history').write_text('another session\n')
        with self.assertRaisesRegex(ValueError, 'cache/history'):
            resume_nextflow_session(self.root)

if __name__ == '__main__': unittest.main()
