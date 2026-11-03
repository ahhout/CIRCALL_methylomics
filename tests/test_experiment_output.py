"""Temporary-only fixtures for experiment ownership and staged publication."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts/lib'))
import experiment_output as adapter
import output_guard as guard

NAME = '2026-11-02_14-00-output-fixture-v1'


class ExperimentOutputs(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.source_mock = patch.object(guard, "require_clean_source", return_value="fixture-commit")
        self.source_mock.start()
        self.addCleanup(self.source_mock.stop)
        self.cli = self.root / "mock_guard_cli.py"
        self.cli.write_text('import sys\nsys.path.insert(0, \'/path/to/hpc-home/analysis/CIRCALL_methylomics/scripts/lib\')\nimport experiment_output, output_guard\n# Explicit test mock: production CLI has no bypass.\noutput_guard.require_clean_source = lambda path: "fixture-commit"\nexperiment_output.main()\n')

    def test_dirty_project_or_experiment_rejects_before_reservation(self):
        for mode in ('project', 'experiment'):
            output = self.root / (NAME + '-' + mode)
            with patch.object(guard, 'require_clean_source', side_effect=(
                    [ValueError('dirty project')] if mode == 'project' else
                    ['project-commit', ValueError('dirty experiment')])):
                with self.assertRaisesRegex(ValueError, 'dirty'):
                    adapter.reserve_bundle(output)
            self.assertFalse(output.exists())
        resource = self.root / 'reference-resource'
        with patch.object(guard, 'require_clean_source', side_effect=ValueError('dirty')):
            with self.assertRaisesRegex(ValueError, 'dirty'):
                adapter.reserve_resource(resource, __file__)
        self.assertFalse(resource.exists())

    def test_existing_objects_and_aliases_never_adopted(self):
        for index, kind in enumerate(['empty', 'nonempty', 'file', 'symlink', 'dangling']):
            output = self.root / (NAME + '-' + str(index))
            if kind in {'empty', 'nonempty'}:
                output.mkdir()
                if kind == 'nonempty': (output / 'original').write_text('unchanged')
            elif kind == 'file': output.write_text('unchanged')
            else: output.symlink_to(self.root if kind == 'symlink' else self.root / 'absent')
            with self.assertRaises(ValueError): adapter.reserve_bundle(output)
            self.assertTrue(os.path.lexists(output))
        alias = self.root / 'alias'; alias.symlink_to(self.root, target_is_directory=True)
        with self.assertRaises(ValueError): adapter.reserve_bundle(alias / NAME)
        with self.assertRaises(ValueError): adapter.reserve_bundle(self.root / 'absent' / '..' / NAME)

    def test_failure_retains_owned_bundle_and_completed_run_rejects_sections(self):
        output = self.root / NAME
        @adapter.owned_execution
        def failed():
            adapter.reserve_bundle(output)
            (output / 'evidence').write_text('retained')
            raise RuntimeError('fixture')
        with self.assertRaises(RuntimeError): failed()
        self.assertEqual(guard.read_record(output)['status'], 'failed')
        self.assertEqual((output / 'evidence').read_text(), 'retained')
        with self.assertRaises(ValueError): adapter.reserve_bundle(output)
        fresh = self.root / (NAME + '-new')
        @adapter.owned_execution
        def succeeded(): adapter.reserve_bundle(fresh)
        succeeded()
        self.assertEqual(guard.read_record(fresh)['status'], 'completed')
        with self.assertRaises(ValueError): adapter.reserve_bundle(fresh / 'new-section', run_root=fresh)

    def test_internal_section_requires_owned_parent_and_cannot_escape(self):
        parent = adapter.reserve_bundle(self.root / NAME)
        child = adapter.reserve_bundle(parent / 'figures', run_root=parent)
        adapter.finish_bundle(child, True)
        self.assertEqual(guard.read_record(child)['status'], 'completed')
        with self.assertRaises(ValueError): adapter.reserve_bundle(self.root / 'escape', run_root=parent)
        unowned = self.root / (NAME + '-unowned'); unowned.mkdir()
        with self.assertRaises(FileNotFoundError): adapter.reserve_bundle(unowned / 'figures', run_root=unowned)

    def test_staging_owner_and_publication_preserve_both_receipts(self):
        owner = adapter.reserve_bundle(self.root / NAME)
        staging = self.root / 'work' / ('.' + NAME + '.running-123') / 'result'
        with patch.dict(os.environ, EXPERIMENT_OWNER_OUTPUT=str(owner)):
            adapter.reserve_bundle(staging, run_name=NAME)
        (staging / 'table.tsv').write_text('fixture')
        adapter.finish_bundle(staging, True)
        result = subprocess.run([sys.executable, str(self.cli), 'publish', '--output', str(owner),
                                 '--staging', str(staging)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((owner / 'table.tsv').read_text(), 'fixture')
        self.assertTrue((owner / 'staging-ownership.json').is_file())
        self.assertEqual(guard.read_record(owner)['status'], 'running')

    def test_prepared_run_checks_snapshots_inputs_and_completion(self):
        output = adapter.reserve_bundle(self.root / NAME)
        (output / 'submitted').mkdir()
        code = output / 'submitted/source.py'; code.write_text('# fixture\n')
        (output / 'source_hashes.json').write_text(json.dumps({'source.py': hashlib.sha256(code.read_bytes()).hexdigest()}))
        (output / 'inputs.json').write_text('{}')
        adapter.freeze_prepared_run(output)
        adapter.require_owned_run(output, source=code)
        unknown = output / "unrecorded.py"; unknown.write_text("# unrecorded")
        with self.assertRaisesRegex(ValueError, "verified submitted snapshot"):
            adapter.require_owned_run(output, source=unknown)
        with patch.object(guard, "require_clean_source", side_effect=["project-commit", "project-commit", ValueError("dirty current helper")]):
            with self.assertRaisesRegex(ValueError, "dirty current helper"):
                adapter.require_owned_run(output, source=__file__)
        original = code.read_text(); code.write_text('# altered\n')
        with self.assertRaises(ValueError): adapter.require_owned_run(output)
        code.write_text(original)
        (output / 'inputs.json').write_text('{"altered": true}')
        with self.assertRaises(ValueError): adapter.require_owned_run(output)
        (output / 'inputs.json').write_text('{}')
        guard.finish(output, True)
        with self.assertRaises(ValueError): adapter.require_owned_run(output)

    def test_completed_ancestor_blocks_staged_bundle_even_with_running_other_owner(self):
        owner = adapter.reserve_bundle(self.root / NAME)
        completed = adapter.reserve_bundle(self.root / (NAME + '-completed'))
        stage = completed / ('.' + NAME + '.running-1')
        stage.mkdir()
        guard.finish(completed, True)
        with patch.dict(os.environ, EXPERIMENT_OWNER_OUTPUT=str(owner)):
            with self.assertRaises(ValueError):
                adapter.reserve_bundle(stage / 'result', run_name=NAME)
        self.assertFalse((stage / 'result').exists())

    def test_report_cannot_write_below_completed_section(self):
        owner = adapter.reserve_bundle(self.root / NAME)
        section = adapter.reserve_bundle(owner / 'figures', run_root=owner)
        adapter.finish_bundle(section, True)
        with patch.dict(os.environ, EXPERIMENT_OWNER_OUTPUT=str(owner)):
            with self.assertRaises(ValueError):
                adapter.prepare_file_outputs([section / 'new.tsv'])
        self.assertFalse((section / 'new.tsv').exists())

    def test_section_completion_cannot_be_reversed(self):
        owner = adapter.reserve_bundle(self.root / NAME)
        section = adapter.reserve_bundle(owner / 'figures', run_root=owner)
        adapter.finish_bundle(section, True)
        original = (section / '.run-ownership.json').read_bytes()
        with self.assertRaises(ValueError): adapter.finish_bundle(section, False)
        self.assertEqual((section / '.run-ownership.json').read_bytes(), original)

    def test_publish_cannot_steal_completed_final_bundle(self):
        owner = adapter.reserve_bundle(self.root / NAME)
        completed = adapter.reserve_bundle(self.root / (NAME + '-completed'))
        (completed / 'payload').write_text('unchanged')
        guard.finish(completed, True)
        before = {p.name: p.read_bytes() for p in completed.iterdir()}
        result = subprocess.run([sys.executable, str(self.cli), 'publish', '--output', str(owner),
                                 '--staging', str(completed)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual({p.name: p.read_bytes() for p in completed.iterdir()}, before)
        self.assertFalse((owner / 'payload').exists())

    def test_nested_completion_finishes_children_before_parent(self):
        owner = self.root / NAME
        @adapter.owned_execution
        def inner():
            adapter.reserve_bundle(owner / 'figures', run_root=owner)
        @adapter.owned_execution
        def outer():
            adapter.reserve_bundle(owner)
            inner()
        outer()
        for path in [owner, owner / 'figures']:
            record = json.loads((path / '.run-ownership.json').read_text())
            self.assertEqual(record['status'], 'completed')

    def test_internal_file_claim_is_exclusive_and_blocks_completed_parent(self):
        owner = adapter.reserve_bundle(self.root / NAME)
        output = owner / 'execution.json'
        command = [sys.executable, str(self.cli), 'claim-file', '--output', str(output)]
        self.assertEqual(subprocess.run(command, capture_output=True).returncode, 0)
        output.write_text('retained')
        self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
        self.assertEqual(output.read_text(), 'retained')
        guard.finish(owner, True)
        command[-1] = str(owner / 'other.json')
        self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
        self.assertFalse((owner / 'other.json').exists())


if __name__ == '__main__':
    unittest.main()
