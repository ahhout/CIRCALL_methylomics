"""Exercise launcher orchestration using tiny fixtures; never invoke conversion."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflows/pat'))
import launch_pat as launcher


class LauncherTests(unittest.TestCase):
    def setUp(self):
        p = patch.object(launcher.output_guard, 'require_clean_source', return_value='a' * 40)
        p.start(); self.addCleanup(p.stop)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'workflows/pat'; self.source.mkdir(parents=True)
        self.results = self.root / 'published'; self.results.mkdir()
        for name in ('prepare_pat.py', 'validate_pat_manifest.py', 'launch_pat.py', 'run_prepare_pat.sh'):
            (self.source / name).write_text('fixture')
        for name in ('environment.yml', 'config/nextflow.config', 'config/params.emseq.integrated.yml',
                     'samplesheets/samplesheet.csv', 'references/emseq_hg38_lambda_puc19/manifest.tsv'):
            p = self.root / name; p.parent.mkdir(parents=True, exist_ok=True); p.write_text('fixture')
        self.input = self.root / 'input'; (self.input / 'bismark/deduplicated').mkdir(parents=True)
        self.metadata = self.root / 'metadata.tsv'; self.metadata.write_text('sample\tpatient\nS1\tP1\n')
        self.sheet = self.root / 'samples.csv'; self.sheet.write_text('sample\nS1\n')
        self.args = argparse.Namespace(results_dir=self.input, metadata=self.metadata, samplesheet=self.sheet,
                                       run_name='2026-09-29_14-35-test', conversion_mode='complete', reuse_dir=[], reuse_sample=[])
        self.calls = []
        for target, value in (('REPO', self.root), ('SOURCE', self.source), ('RESULTS', self.results)):
            p = patch.object(launcher, target, value); p.start(); self.addCleanup(p.stop)
        p = patch.dict(os.environ, {'SLURM_JOB_ID': '123', 'SLURM_CPUS_PER_TASK': '8'}); p.start(); self.addCleanup(p.stop)
        p = patch.object(launcher.shutil, 'which', return_value='/fixture/tool'); p.start(); self.addCleanup(p.stop)
        p = patch.object(launcher.subprocess, 'check_output', return_value='fixture-revision\n'); p.start(); self.addCleanup(p.stop)

    def fake_command(self, argv, log=None):
        argv = [str(x) for x in argv]; self.calls.append(argv)
        if '--output-dir' in argv:
            bundle = Path(argv[argv.index('--output-dir') + 1]); bundle.mkdir()
            (bundle / 'PAT_COMPLETED').write_text('fixture')
        elif argv[0] == 'mv':
            Path(argv[-2]).rename(argv[-1])
        else:
            self.assertTrue(Path(argv[1]).is_file())
            self.assertTrue(Path(argv[argv.index('--samplesheet') + 1]).is_file())
            self.assertTrue(Path(argv[argv.index('--pat-dir') + 1]).is_dir())
            log.write_text('fixture PASS\n')

    def test_required_arguments_and_modes(self):
        with self.assertRaises(SystemExit): launcher.parse_args([])
        with self.assertRaises(SystemExit):
            launcher.parse_args(['--results-dir','x','--metadata','x','--samplesheet','x',
                                 '--run-name','x','--conversion-mode','rebuild'])

    def test_timestamp_required_before_creating_work(self):
        for name in ('2026-09-29-test', '2026-02-30_10-00-test',
                     '2026-09-29_24-00-test', '2026-09-29_10-60-test'):
            with self.subTest(name=name):
                self.args.run_name = name
                with self.assertRaises(ValueError): launcher.launch(self.args)
                self.assertEqual(list(self.results.iterdir()), [])

    def test_complete_generates_and_validates_before_and_after_publication(self):
        with patch.object(launcher, 'command', side_effect=self.fake_command): launcher.launch(self.args)
        first = self.calls[0]
        self.assertEqual(first[first.index('--conversion-mode') + 1], 'complete')
        self.assertEqual(first[first.index('--threads') + 1], '8')
        self.assertNotIn('--reuse-dir', first)
        final = self.results / '2026-09-29_14-35-test-pat-v1'
        self.assertTrue((final / 'submitted.sha256').is_file())
        self.assertTrue((final / 'validation_after_publish.log').is_file())
        self.assertEqual(len(self.calls), 3)
        timing=json.loads((final/'publication_execution.json').read_text())
        self.assertEqual(timing['status'], 'completed')
        self.assertIn('started_utc', timing);self.assertIn('completed_utc', timing)

    def test_reuse_and_explicit_partial_reuse_forwarding(self):
        bundle = self.root / 'reuse'; bundle.mkdir()
        for mode in ('reuse', 'complete'):
            self.args.run_name = '2026-09-29_14-35-' + mode; self.args.conversion_mode = mode
            self.args.reuse_dir = [bundle]
            self.args.reuse_sample = ['S1'] if mode == 'complete' else []
            self.calls.clear()
            with patch.object(launcher, 'command', side_effect=self.fake_command): launcher.launch(self.args)
            first = self.calls[0]
            self.assertEqual(first[first.index('--reuse-dir') + 1], str(bundle))
            self.assertEqual('--reuse-sample' in first, mode == 'complete')

    def test_mismatched_samples_fail_before_work(self):
        self.sheet.write_text('sample\nS2\n')
        with self.assertRaisesRegex(ValueError, 'identical samples'): launcher.launch(self.args)
        self.assertEqual(list(self.results.iterdir()), [])

    def test_existing_output_and_broken_symlink_refused(self):
        final = self.results / '2026-09-29_14-35-test-pat-v1'; final.symlink_to(self.root / 'missing')
        with self.assertRaisesRegex(ValueError, 'destination exists'): launcher.launch(self.args)
        self.assertTrue(final.is_symlink())

    def test_existing_work_refused(self):
        (self.results / '.2026-09-29_14-35-test-pat-v1.running-123').mkdir()
        with self.assertRaisesRegex(ValueError, 'destination exists'): launcher.launch(self.args)

    def test_concurrent_lock_refused(self):
        with (self.results / '2026-09-29_14-35-test-pat-v1.lock').open('a') as handle:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaises(BlockingIOError): launcher.launch(self.args)
        self.assertFalse((self.results / '.2026-09-29_14-35-test-pat-v1.running-123').exists())

    def test_failure_propagates_and_retains_workspace_without_publication(self):
        with patch.object(launcher, 'command', side_effect=subprocess.CalledProcessError(7, 'fixture')):
            with self.assertRaises(subprocess.CalledProcessError): launcher.launch(self.args)
        self.assertTrue((self.results / '.2026-09-29_14-35-test-pat-v1.running-123/submitted').is_dir())
        timing=json.loads((self.results/'.2026-09-29_14-35-test-pat-v1.running-123/publication_execution.json').read_text())
        self.assertIn('started_utc',timing);self.assertIn('ended_utc',timing)
        self.assertNotIn('completed_utc',timing)
        self.assertEqual(launcher.output_guard.read_record(self.results / '2026-09-29_14-35-test-pat-v1')['status'], 'failed')

    def test_late_collision_preserves_destination_and_retains_work(self):
        final=self.results/'2026-09-29_14-35-test-pat-v1'
        def late_collision(argv, log=None):
            self.fake_command(argv,log)
            if log and log.name == 'validation_before_publish.log':
                (final/'PAT_COMPLETED').write_text('original')
        with patch.object(launcher,'command',side_effect=late_collision):
            with self.assertRaisesRegex(FileExistsError,'conflicts'): launcher.launch(self.args)
        self.assertEqual((final/'PAT_COMPLETED').read_text(),'original')
        work=self.results/'.2026-09-29_14-35-test-pat-v1.running-123/bundle'
        self.assertTrue(work.is_dir())
        timing=json.loads((work/'publication_execution.json').read_text())
        self.assertNotIn('completed_utc',timing)
        self.assertEqual(timing['status'],'failed_or_interrupted')

    def test_post_publication_validation_failure_has_no_success_timestamp(self):
        def fail_after_publication(argv, log=None):
            if log and log.name == 'validation_after_publish.log':
                raise subprocess.CalledProcessError(9, 'validator')
            self.fake_command(argv, log)
        with patch.object(launcher, 'command', side_effect=fail_after_publication):
            with self.assertRaises(subprocess.CalledProcessError): launcher.launch(self.args)
        final=self.results/'2026-09-29_14-35-test-pat-v1'
        timing=json.loads((final/'publication_execution.json').read_text())
        self.assertEqual(timing['status'],'failed_or_interrupted')
        self.assertIn('ended_utc',timing)
        self.assertNotIn('completed_utc',timing)

    def test_validation_failure_does_not_publish(self):
        def fail_validation(argv, log=None):
            if log: raise subprocess.CalledProcessError(8, 'validator')
            self.fake_command(argv, log)
        with patch.object(launcher, 'command', side_effect=fail_validation):
            with self.assertRaises(subprocess.CalledProcessError): launcher.launch(self.args)
        self.assertTrue((self.results / '.2026-09-29_14-35-test-pat-v1.running-123/bundle').is_dir())
        self.assertEqual(launcher.output_guard.read_record(self.results / '2026-09-29_14-35-test-pat-v1')['status'], 'failed')


if __name__ == '__main__': unittest.main()
