import csv
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / 'archive/downstream/dmrichr'
spec = importlib.util.spec_from_file_location('dmrichr_prepare', WORKFLOW / 'prepare_inputs.py')
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class DmrichrInputsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.dirs = [self.root / 'old', self.root / 'new']
        for directory in self.dirs:
            directory.mkdir()
        self.design = self.root / 'sample_info.tsv'
        self.metadata = self.root / 'metadata.tsv'
        with self.design.open('w') as design, self.metadata.open('w') as metadata:
            d = csv.writer(design, delimiter='\t'); m = csv.writer(metadata, delimiter='\t')
            d.writerow(['Name', 'Compartment', 'Patient']); m.writerow(['sample', 'patient'])
            for patient in range(4):
                for group in ['A_PB', 'B_BM']:
                    sample = f'S{patient}{group[0]}'
                    d.writerow([sample, group, f'P{patient}']); m.writerow([sample, f'P{patient}'])
                    with gzip.open(self.dirs[patient == 3] / f'{sample}.CpG_report.txt.gz', 'wt') as f:
                        f.write('chr1\t100\t+\t3\t2\tCG\tCGA\n')
        self.args = SimpleNamespace(sample_info=self.design, metadata=self.metadata,
            reports_dir=self.dirs, work_dir=self.root/'work', analysis_id='paired-four',
            genome='hg38', test_covariate='Compartment', match_covariate='Patient',
            runtime_label='fixture', rng_seed=20260821)

    def test_eight_samples_use_two_explicit_sources_without_copying(self):
        prepare.run(self.args)
        manifest = json.loads((self.args.work_dir/'input_manifest.json').read_text())
        self.assertEqual(len(manifest['reports']), 8)
        self.assertEqual(manifest['reported_direction_expected'], 'B_BM minus A_PB')
        for record in manifest['reports']:
            staged = self.args.work_dir / (record['sample'] + '_integrated.CpG_report.txt.gz')
            self.assertTrue(staged.is_symlink())
            self.assertEqual(str(staged.resolve()), record['source'])
        self.assertLessEqual(manifest['preparation_started_utc'], manifest['preparation_completed_utc'])

    def test_missing_new_cohort_fails_before_staging(self):
        self.args.reports_dir = self.dirs[:1]
        with self.assertRaisesRegex(ValueError, 'exactly one CpG report'):
            prepare.run(self.args)
        self.assertFalse(self.args.work_dir.exists())

    def test_ambiguous_sample_fails_before_staging(self):
        (self.dirs[1]/'S0A.CpG_report.txt.gz').symlink_to(self.dirs[0]/'S0A.CpG_report.txt.gz')
        with self.assertRaisesRegex(ValueError, 'found 2'):
            prepare.run(self.args)
        self.assertFalse(self.args.work_dir.exists())

    def test_duplicate_patient_group_is_rejected(self):
        self.design.write_text(self.design.read_text().replace('S0B\tB_BM', 'S0B\tA_PB'))
        with self.assertRaisesRegex(ValueError, 'exactly one sample per group'):
            prepare.run(self.args)

    def test_execution_records_preserve_original_start_across_reporting(self):
        path = self.root/'execution.json'
        command = [sys.executable, str(WORKFLOW/'record_execution.py'), '--manifest', str(path),
                   '--analysis-id', 'fixture', '--run-name', '2026-09-29_10-00-fixture-v1']
        for phase, job, start in [('core', '1', '2026-09-29T10:00:00Z'),
                                  ('reporting', '2', '2026-09-29T11:00:00Z')]:
            for event in ['started', 'completed']:
                subprocess.run(command + ['--job-id', job, '--phase', phase, '--event', event,
                                          '--started-utc', start], check=True)
        result = json.loads(path.read_text())
        self.assertEqual(result['started_utc'], '2026-09-29T10:00:00Z')
        self.assertEqual(len(result['attempts']), 2)
        self.assertTrue(result['completed_utc'])


if __name__ == '__main__':
    unittest.main()
