import gzip
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import methylseq_smoke as smoke
import test_methylseq_launch as fixtures
launch, NAME = fixtures.launch, fixtures.NAME


class SmokeTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.LaunchTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.project = self.fixture.project
        self.reads = []
        for mate in (1, 2):
            p = self.fixture.base / f'reads_{mate}.fq.gz'
            with gzip.open(p, 'wb') as handle:
                for i in range(4):
                    handle.write(f'@pair{i}/{mate}\nACGT\n+\n!!!!\n'.encode())
            self.reads.append(p)
        (self.project / 'samplesheets/samplesheet.csv').write_text(
            'sample,fastq_1,fastq_2,genome\nfixture,' + ','.join(map(str, self.reads)) + ',\n')
        (self.project / 'config/methylseq_smoke_resources.config').write_text("process.queue = 'shortq'\n")
        self.fixture.commit()

    def prepare(self):
        return smoke.prepare(self.project, NAME, 2)

    def test_prefix_preserves_pair_identity_and_source_bytes(self):
        before = [p.read_bytes() for p in self.reads]
        manifest = self.prepare()
        value = smoke.validate(self.project, manifest)
        self.assertEqual(value['pairs'], 2)
        for mate in (1, 2):
            with gzip.open(manifest.parent / f'reads_R{mate}.fastq.gz', 'rt') as handle:
                lines = handle.readlines()
            self.assertEqual(len(lines), 8)
            self.assertEqual(lines[0].strip(), f'@pair0/{mate}')
        self.assertEqual(before, [p.read_bytes() for p in self.reads])
        with self.assertRaises(FileExistsError): self.prepare()

    def test_changed_manifest_or_payload_rejected(self):
        manifest = self.prepare()
        original = manifest.read_bytes()
        manifest.write_text('{}')
        with self.assertRaisesRegex(ValueError, 'altered'): smoke.validate(self.project, manifest)
        manifest.write_bytes(original)
        (manifest.parent / 'reads_R1.fastq.gz').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'artifact changed'): smoke.validate(self.project, manifest)

    def test_invalid_pairing_retains_failed_reservation(self):
        with gzip.open(self.reads[1], 'wb') as handle:
            handle.write(b'@wrong/2\nACGT\n+\n!!!!\n')
        with self.assertRaisesRegex(ValueError, 'out of sync'): self.prepare()
        bundle = self.project / 'runs/methylseq/smoke_inputs' / NAME
        self.assertEqual(smoke.guard.read_record(bundle)['status'], 'failed')
        with self.assertRaises(FileExistsError): self.prepare()

    def test_smoke_uses_same_reference_new_session_and_owned_outputs(self):
        manifest = self.prepare()
        params = (self.project / 'config/params.emseq.integrated.yml').read_bytes()
        output = launch.execute(self.project, NAME, runner=self.fixture.runner, smoke_inputs=manifest)
        argv = self.fixture.calls[0]
        self.assertNotIn('-resume', argv)
        self.assertIn('smoke_resources.config', launch.option(argv, '-C'))
        self.assertEqual(launch.option(argv, '-work-dir'), str(self.fixture.paths['work_root']))
        submitted = output / 'pipeline_info/submitted'
        self.assertEqual((submitted / 'params.emseq.integrated.yml').read_bytes(), params)
        self.assertEqual((submitted / 'samplesheet.csv').read_bytes(), (manifest.parent / 'samplesheet.csv').read_bytes())
        with self.assertRaises((ValueError, FileExistsError)):
            launch.execute(self.project, NAME, runner=self.fixture.runner, smoke_inputs=manifest)

    def test_mode_is_frozen_and_scheduler_receives_manifest(self):
        manifest = self.prepare()
        _, _, identity = launch.preflight(self.project, NAME, manifest)
        _, command = launch.submission_command(self.project, NAME, False, identity, manifest)
        self.assertIn('CIRCALL_SMOKE_INPUTS=' + str(manifest), command[1])
        self.assertIn('--partition=shortq', command)
        self.fixture.pipeline_status = 1
        with self.assertRaises(ValueError):
            launch.execute(self.project, NAME, runner=self.fixture.runner, smoke_inputs=manifest)
        self.fixture.pipeline_status = 0
        with self.assertRaisesRegex(ValueError, 'identity'):
            launch.execute(self.project, NAME, True, runner=self.fixture.runner)
        launch.execute(self.project, NAME, True, runner=self.fixture.runner, smoke_inputs=manifest)
        self.assertIn('-resume', self.fixture.calls[1])

    def test_bound_and_external_manifest_fail_before_creation(self):
        with self.assertRaises(ValueError): smoke.prepare(self.project, NAME, 10001)
        self.assertFalse((self.project / 'runs/methylseq/smoke_inputs' / NAME).exists())
        manifest = self.prepare()
        external = self.fixture.base / 'manifest.json'; external.write_bytes(manifest.read_bytes())
        with self.assertRaisesRegex(ValueError, 'outside'): smoke.validate(self.project, external)


if __name__ == '__main__':
    unittest.main()
