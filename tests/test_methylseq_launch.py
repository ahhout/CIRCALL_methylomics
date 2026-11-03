#!/usr/bin/env python3
"""No Nextflow/Slurm jobs: all execution, sessions and input data are fixtures."""
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import uuid

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('methylseq_launch', ROOT / 'scripts/methylseq_launch.py')
launch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launch)
NAME = '2026-11-02_12-30-emseq-integrated-fixture-v1'


class LaunchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='methylseq-launch-fixture-')
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.project = self.base / 'project'
        self.pipeline = self.base / 'pipeline'
        self.project.mkdir()
        self.pipeline.mkdir()
        (self.pipeline / 'main.nf').write_text('// fixture pipeline\n')
        (self.pipeline / 'nextflow.config').write_text('// fixture config\n')
        pipeline_commit = self.init_git(self.pipeline)
        paths = {
            'schema': 1, 'pipeline_root': str(self.pipeline),
            'runtime_root': str(self.project / 'runs/methylseq'),
            'work_root': str(self.base / 'work'), 'maintenance_lock': str(self.base / 'shared.lock'),
            'results_parent': str(self.base / 'results'),
            'reference_root': str(self.project / 'references/emseq_hg38_lambda_puc19'),
            'container_cache': str(self.base / 'containers')}
        self.paths = {key: Path(value) for key, value in paths.items() if key != 'schema'}
        for key in ('work_root', 'results_parent', 'container_cache', 'reference_root'):
            self.paths[key].mkdir(parents=True)
        for sub in ('config', 'scripts', 'samplesheets', 'outputs'):
            (self.project / sub).mkdir()
        reference = self.paths['reference_root']
        (reference / 'manifest.tsv').write_text('fixture\tidentity\n')
        for conversion in ('CT_conversion', 'GA_conversion'):
            target = reference / 'BismarkIndex/Bisulfite_Genome' / conversion
            target.mkdir(parents=True)
            (target / 'fixture.1.bt2').write_text('fixture index')
        (reference / 'BismarkIndex/genome.fa').write_text('>phage_lambda\nAC\n>plasmid_puc19c\nGT\n')
        (self.project / 'config/methylseq_paths.json').write_text(json.dumps(paths))
        (self.project / 'config/methylseq_pipeline.json').write_text(json.dumps(dict(
            schema=1, pipeline='nf-core/methylseq', release='4.2.0', commit=pipeline_commit,
            upstream_release_commit=pipeline_commit, fork_baseline_commit=pipeline_commit)))
        (self.project / 'config/nextflow.config').write_text("process.executor = 'slurm'\n")
        (self.project / 'config/params.emseq.integrated.yml').write_text(
            f"input: {self.project / 'samplesheets/samplesheet.csv'}\n"
            f"fasta: {reference / 'BismarkIndex/genome.fa'}\n"
            f"bismark_index: {reference / 'BismarkIndex'}\naligner: bismark\nem_seq: true\n")
        (self.project / 'run_methylseq_integrated.sh').write_text('#!/bin/bash\nset -euo pipefail\n')
        fastq = self.base / 'fixture.fastq'
        fastq.write_text('@read\nAC\n+\n!!\n')
        (self.project / 'samplesheets/samplesheet.csv').write_text(
            f'sample,fastq_1,fastq_2,genome\nfixture,{fastq},,\n')
        (self.project / '.gitignore').write_text('/runs/\n/outputs/\n')
        (self.project / 'scripts/summarize_emseq_controls.py').write_text('# fixture\n')
        self.init_git(self.project)
        self.calls = []
        self.sessions = {}
        self.pipeline_status = 0
        self.addCleanup(patch.stopall)
        patch.object(launch, 'check_versions', return_value={'nextflow': '25.10.2', 'apptainer': '1.3.6'}).start()
        patch.dict(os.environ, {'SLURM_JOB_ID': '123456', 'NXF_CONFIG': '', 'NXF_PARAMS': '', 'NXF_CLI': ''}).start()

    def init_git(self, repo):
        def command(*args):
            return subprocess.check_output(['git', '-C', str(repo), *args], stderr=subprocess.DEVNULL, text=True).strip()
        command('init', '-q')
        command('config', 'user.name', 'Fixture')
        command('config', 'user.email', 'fixture@example.invalid')
        command('add', '.')
        command('commit', '-qm', 'Fixture source')
        return command('rev-parse', 'HEAD')

    def commit(self):
        subprocess.run(['git', '-C', str(self.project), 'add', '.'], check=True)
        subprocess.run(['git', '-C', str(self.project), 'commit', '-qm', 'Fixture change'], check=True)

    def runner(self, command, **kwargs):
        self.calls.append(command)
        if command[0] == 'nextflow':
            output = launch.option(command, '--outdir')
            runtime = Path(kwargs['cwd'])
            session = self.sessions.setdefault(output, str(uuid.uuid4()))
            cache = runtime / '.nextflow/cache' / session / 'db'
            cache.mkdir(parents=True, exist_ok=True)
            (cache / 'CURRENT').write_text('fixture cache metadata')
            with (runtime / '.nextflow/history').open('a') as history:
                history.write('start\tduration\trun\tOK\t1\t' + session + '\t' + shlex.join(command) + '\n')
            return subprocess.CompletedProcess(command, self.pipeline_status)
        control = Path(launch.option(command, '--output'))
        control.parent.mkdir(parents=True, exist_ok=True)
        control.write_text('fixture control QC\n')
        return subprocess.CompletedProcess(command, 0)

    def run_fixture(self, resume=False, name=NAME):
        return launch.execute(self.project, name, resume, runner=self.runner)

    def test_new_session_no_resume_and_completed_is_immutable(self):
        output = self.run_fixture()
        self.assertNotIn('-resume', self.calls[0])
        self.assertIn('-C', self.calls[0])
        self.assertEqual(launch.option(self.calls[0], '--custom_config_base'), 'false')
        self.assertEqual(launch.option(self.calls[0], '-profile'), 'apptainer')
        self.assertEqual(launch.option(self.calls[0], '-work-dir'), str(self.paths['work_root']))
        self.assertEqual(Path(launch.option(self.calls[0], '--input')).parent, output / 'pipeline_info/submitted')
        record = launch.output_guard.read_record(output)
        self.assertEqual(record['status'], 'completed')
        before = (output / '.run-ownership.json').read_bytes()
        with self.assertRaises((ValueError, FileExistsError)):
            self.run_fixture()
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            self.run_fixture(resume=True)
        self.assertEqual(before, (output / '.run-ownership.json').read_bytes())

    def test_failed_exact_resume_keeps_initial_snapshot_and_timestamp(self):
        self.pipeline_status = 1
        with self.assertRaisesRegex(ValueError, 'status 1'):
            self.run_fixture()
        output = self.paths['results_parent'] / NAME
        before = launch.output_guard.read_record(output)
        snapshot = output / 'pipeline_info/submitted/run_manifest.json'
        original = snapshot.read_bytes()
        self.pipeline_status = 0
        self.run_fixture(resume=True)
        after = launch.output_guard.read_record(output)
        self.assertEqual(before['started_utc'], after['started_utc'])
        self.assertEqual(before['session_id'], after['session_id'])
        self.assertEqual(original, snapshot.read_bytes())
        self.assertEqual(len(after['attempts']), 2)
        self.assertEqual(launch.option(self.calls[1], '-resume'), before['session_id'])
        self.assertEqual(len(list((output / 'pipeline_info/submitted/attempts').iterdir())), 2)

    def fail_run(self):
        self.pipeline_status = 1
        with self.assertRaises(ValueError):
            self.run_fixture()
        return self.paths['results_parent'] / NAME

    def test_altered_snapshot_is_rejected_before_execution(self):
        output = self.fail_run()
        (output / 'pipeline_info/submitted/samplesheet.csv').write_text('changed')
        before = (output / '.run-ownership.json').read_bytes()
        with self.assertRaisesRegex(ValueError, 'snapshot changed'):
            self.run_fixture(resume=True)
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(before, (output / '.run-ownership.json').read_bytes())

    def test_missing_and_conflicting_sessions_rejected(self):
        output = self.fail_run()
        history = self.paths['runtime_root'] / '.nextflow/history'
        original = history.read_text()
        history.write_text(original + original.replace(self.sessions[str(output)], str(uuid.uuid4())))
        with self.assertRaisesRegex(ValueError, 'conflicting session'):
            self.run_fixture(resume=True)
        history.write_text(original)
        session = self.sessions[str(output)]
        (self.paths['runtime_root'] / '.nextflow/cache' / session / 'db/CURRENT').unlink()
        with self.assertRaisesRegex(ValueError, 'metadata is missing'):
            self.run_fixture(resume=True)

    def test_original_manifest_and_checksum_inventory_are_immutable(self):
        output = self.fail_run()
        submitted = output / 'pipeline_info/submitted'
        for name in ('SHA256SUMS', 'run_manifest.json'):
            path = submitted / name
            original = path.read_bytes()
            path.write_text('altered')
            with self.assertRaisesRegex(ValueError, 'changed'):
                self.run_fixture(resume=True)
            path.write_bytes(original)
        self.assertEqual(len(self.calls), 1)

    def test_preflight_destination_refuses_existing_empty_directory_without_mutation(self):
        paths, files, identity = launch.preflight(self.project, NAME)
        output = self.paths['results_parent'] / NAME
        output.mkdir()
        with self.assertRaisesRegex(ValueError, 'already exists'):
            launch.check_destination(paths, files, identity, NAME, False)
        self.assertEqual(list(output.iterdir()), [])

    def test_changed_work_in_history_rejected(self):
        self.fail_run()
        history = self.paths['runtime_root'] / '.nextflow/history'
        history.write_text(history.read_text().replace(str(self.paths['work_root']), str(self.base / 'different')))
        with self.assertRaisesRegex(ValueError, 'work-root mismatch'):
            self.run_fixture(resume=True)

    def test_dirty_source_and_wrong_pipeline_pin_rejected(self):
        (self.project / 'untracked-source.py').write_text('changed')
        with self.assertRaisesRegex(ValueError, 'Commit reviewed'):
            launch.preflight(self.project, NAME)
        (self.project / 'untracked-source.py').unlink()
        pin = self.project / 'config/methylseq_pipeline.json'
        record = json.loads(pin.read_text())
        record['commit'] = '0' * 40
        pin.write_text(json.dumps(record))
        self.commit()
        with self.assertRaisesRegex(ValueError, 'pinned commit'):
            launch.preflight(self.project, NAME)

    def test_output_and_external_config_overrides_rejected(self):
        params = self.project / 'config/params.emseq.integrated.yml'
        original = params.read_text()
        for key in ('outdir', 'custom_config_base', 'trace_report_suffix'):
            params.write_text(original + f'{key}: arbitrary\n')
            with self.assertRaisesRegex(ValueError, 'Conflicting output'):
                launch.preflight(self.project, NAME)
        params.write_text(original)
        (self.project / 'config/nextflow.config').write_text("process.publishDir = '/outside'\n")
        with self.assertRaisesRegex(ValueError, 'override publication'):
            launch.preflight(self.project, NAME)
        self.assertFalse((self.paths['results_parent'] / NAME).exists())

    def test_lock_contention_blocks_before_reservation(self):
        with launch.execution_lock(self.paths['maintenance_lock']):
            with self.assertRaisesRegex(ValueError, 'in use'):
                self.run_fixture()
        self.assertFalse((self.paths['results_parent'] / NAME).exists())
        self.assertEqual(self.calls, [])

    def test_conflicting_navigation_link_preserved_and_failure_retained(self):
        target = self.project / 'outputs' / NAME
        target.symlink_to(self.base / 'somewhere-else')
        with self.assertRaisesRegex(ValueError, 'another output'):
            self.run_fixture()
        self.assertEqual(target.readlink(), self.base / 'somewhere-else')
        record = launch.output_guard.read_record(self.paths['results_parent'] / NAME)
        self.assertEqual(record['status'], 'failed')
        self.assertIsNone(record['session_id'])
        self.assertEqual(self.calls, [])
        with self.assertRaisesRegex(ValueError, 'lacks a recorded session'):
            self.run_fixture(resume=True)

    def test_existing_empty_file_and_dangling_links_rejected(self):
        for kind in ('empty', 'file', 'symlink'):
            name = NAME.replace('fixture', kind)
            output = self.paths['results_parent'] / name
            if kind == 'empty':
                output.mkdir()
            elif kind == 'file':
                output.write_text('keep')
            else:
                output.symlink_to(self.base / 'missing')
            with self.assertRaises(FileExistsError):
                self.run_fixture(name=name)
            self.assertTrue(os.path.lexists(output))
        self.assertEqual(self.calls, [])

    def test_distinct_runs_separate_sessions_and_outputs(self):
        first = self.run_fixture()
        second = self.run_fixture(name=NAME.replace('fixture', 'another'))
        self.assertNotEqual(first, second)
        self.assertNotEqual(launch.output_guard.read_record(first)['session_id'],
                            launch.output_guard.read_record(second)['session_id'])
        self.assertNotIn('-resume', self.calls[2])

    def test_submission_records_commits_with_logs_outside_output(self):
        paths, files, identity = launch.preflight(self.project, NAME)
        logs, command = launch.submission_command(self.project, NAME, False, identity)
        self.assertTrue(logs.is_relative_to(self.project / 'runs/methylseq/scheduler'))
        self.assertFalse(logs.exists())
        self.assertEqual(command[0], 'sbatch')
        self.assertIn('CIRCALL_PROJECT_COMMIT=' + identity['project_commit'], command[1])
        self.assertIn('CIRCALL_PIPELINE_COMMIT=' + identity['pipeline_commit'], command[1])
        self.assertTrue(any('%j.out' in item for item in command))
        self.assertNotIn(str(paths['results_parent']), ' '.join(command))

    def test_invalid_calendar_and_noncanonical_parent_rejected(self):
        for name in (NAME.replace('11-02', '02-30'), NAME.replace('12-30', '25-30'), '../' + NAME):
            with self.assertRaises(ValueError):
                launch.preflight(self.project, name)
        alias = self.base / 'alias'
        alias.symlink_to(self.paths['results_parent'])
        config = self.project / 'config/methylseq_paths.json'
        paths = json.loads(config.read_text())
        paths['results_parent'] = str(alias)
        config.write_text(json.dumps(paths))
        with self.assertRaisesRegex(ValueError, 'noncanonical'):
            launch.preflight(self.project, NAME)


if __name__ == '__main__':
    unittest.main()
