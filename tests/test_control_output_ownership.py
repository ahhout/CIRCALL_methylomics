"""Direct control-QC CLI fixtures cannot overwrite completed or unowned results."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import unittest
import test_methylseq_launch as fixtures


class ControlOutputTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.LaunchTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.project = self.fixture.project
        shutil.copy2(fixtures.ROOT / 'scripts/summarize_emseq_controls.py', self.project / 'scripts/summarize_emseq_controls.py')
        (self.project / 'scripts/lib').mkdir()
        shutil.copy2(fixtures.ROOT / 'scripts/lib/output_guard.py', self.project / 'scripts/lib/output_guard.py')
        self.fixture.commit()
        self.output = self.fixture.fail_run()
        self.record = fixtures.launch.output_guard.read_record(self.output)
        self.attempt = self.record['attempts'][-1]['attempt_id']
        self.report = self.output / 'control_qc/control_qc.tsv'
        coverage = self.output / 'bismark/methylation_calls/methylation_coverage'
        coverage.mkdir(parents=True)
        with gzip.open(coverage / 'fixture.coverage.cov.gz', 'wt') as handle:
            handle.write('phage_lambda\t1\t1\t0\t0\t100\nplasmid_puc19c\t1\t1\t100\t100\t0\n')
        self.command = [sys.executable, str(self.project / 'scripts/summarize_emseq_controls.py'),
                        '--coverage-dir', str(coverage), '--samplesheet', str(self.output / 'pipeline_info/submitted/samplesheet.csv'),
                        '--output', str(self.report), '--owned-output-root', str(self.output), '--attempt-id', self.attempt]

    def invoke(self):
        return subprocess.run(self.command, env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'),
                              capture_output=True, text=True, timeout=20)

    def test_direct_completed_output_is_immutable(self):
        self.report.parent.mkdir()
        self.report.write_text('prior result')
        fixtures.launch.output_guard.finish(self.output, True)
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('immutable', result.stderr)
        self.assertEqual(self.report.read_text(), 'prior result')

    def test_unowned_historical_directory_is_rejected(self):
        (self.output / '.run-ownership.json').unlink()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.report.exists())

    def test_symlink_output_is_rejected(self):
        fixtures.launch.output_guard.update(self.output, status='running')
        target = self.fixture.base / 'untouched.tsv'
        target.write_text('keep')
        self.report.parent.mkdir()
        self.report.symlink_to(target)
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(target.read_text(), 'keep')

    def test_owned_running_qc_succeeds_and_records_timestamps(self):
        fixtures.launch.output_guard.update(self.output, status='running')
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('100.000000', self.report.read_text())
        provenance = json.loads((self.output / 'pipeline_info/submitted/attempts' / self.attempt / 'control-qc.json').read_text())
        self.assertIn('started_utc', provenance)
        self.assertIn('completed_utc', provenance)
        second = self.invoke()
        self.assertNotEqual(second.returncode, 0)


if __name__ == '__main__':
    unittest.main()
