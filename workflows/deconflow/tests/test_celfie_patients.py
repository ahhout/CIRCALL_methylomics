"""Patient count integrity, native formatting, observer and resume boundaries."""
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import ATLAS_SHA, read_json, sha, write_json
from tools import celfie_patients as route
from tools.celfie_diagnostics import observe


def counts(path, labels, records=None):
    header = ['chr', 'start', 'end'] + [f'{name}_{kind}' for name in labels for kind in ('METH', 'DEPTH')]
    records = records or [['chr1', '10', '20', *['0', '0'] * len(labels)],
                          ['chr2', '30', '40', *['3', '5'] * len(labels)]]
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(records)


class PatientCounts(unittest.TestCase):
    def test_zero_depth_is_preserved_and_native_layout_retains_exact_pairs(self):
        coordinates = [('chr1', '10', '20'), ('chr2', '30', '40')]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts(root / 'reference.tsv', ['A', 'B'])
            counts(root / 'patients.tsv', ['sample_two', 'sample_one'])
            reference, _ = route.count_table(root / 'reference.tsv', coordinates, ['A', 'B'])
            patient, _ = route.count_table(root / 'patients.tsv', coordinates)
            text = route.native_table(reference, patient, ['A', 'B'], ['sample_one'])
            lines = text.splitlines()
            self.assertEqual(lines[0].split('\t'), ['CHR', 'START', 'END', 'sample_one_meth',
                             'sample_one_depth', 'REF_CHR', 'REF_START', 'REF_END',
                             'A_meth', 'A_depth', 'B_meth', 'B_depth'])
            self.assertEqual(lines[1].split('\t'), ['chr1', '10', '20', '0', '0',
                             'chr1', '10', '20', '0', '0', '0', '0'])
            self.assertEqual(len(lines), 3)

    def test_rejects_counts_coordinate_and_reference_label_drift(self):
        coordinates = [('chr1', '10', '20'), ('chr2', '30', '40')]
        for change in ('negative', 'fractional', 'missing', 'meth_gt_depth', 'coordinate', 'order'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'counts.tsv'
                records = [['chr1', '10', '20', '0', '0'], ['chr2', '30', '40', '3', '5']]
                if change == 'negative':
                    records[0][3] = '-1'
                elif change == 'fractional':
                    records[0][3] = '0.0'
                elif change == 'missing':
                    records[0][3] = ''
                elif change == 'meth_gt_depth':
                    records[0][3] = '1'
                elif change == 'coordinate':
                    records[0][1] = '11'
                counts(path, ['A'], records)
                with self.assertRaises(ValueError):
                    route.count_table(path, coordinates, ['B'] if change == 'order' else ['A'])

    def test_gate_binds_actual_reference_and_optional_patient_count_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts(root / 'reference.tsv', ['A'])
            counts(root / 'patients.tsv', ['sample'])
            gate = dict(status='passed', atlas_sha256=ATLAS_SHA,
                        reference_counts_sha256=sha(root / 'reference.tsv'),
                        patient_counts_sha256=sha(root / 'patients.tsv'))
            write_json(root / 'gate.json', gate)
            route.gate_reference(root / 'gate.json', root / 'reference.tsv', root / 'patients.tsv')
            with (root / 'reference.tsv').open('a') as handle:
                handle.write('\n')
            with self.assertRaisesRegex(ValueError, 'Reference counts differ'):
                route.gate_reference(root / 'gate.json', root / 'reference.tsv', root / 'patients.tsv')

    def test_changed_run_local_configuration_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / 'celfie.config'
            route.immutable_text(config, "executor = 'slurm'\n")
            config.write_text("executor = 'local'\n")
            with self.assertRaisesRegex(ValueError, 'changed run-local'):
                route.immutable_text(config, "executor = 'slurm'\n")

    def test_run_lock_rejects_a_simultaneous_attempt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with route.run_lock(root):
                with self.assertRaisesRegex(ValueError, 'already running'):
                    with route.run_lock(root):
                        self.fail('Second attempt acquired the same run lock')
            with route.run_lock(root):
                pass

    def test_sample_with_no_observed_counts_is_rejected_before_fitting(self):
        coordinates = [('chr1', '10', '20'), ('chr2', '30', '40')]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts(root / 'reference.tsv', ['A', 'B'])
            counts(root / 'patients.tsv', ['sample'], [['chr1', '10', '20', '0', '0'],
                                                      ['chr2', '30', '40', '0', '0']])
            reference, _ = route.count_table(root / 'reference.tsv', coordinates, ['A', 'B'])
            patients, _ = route.count_table(root / 'patients.tsv', coordinates)
            with self.assertRaisesRegex(ValueError, 'No observed counts'):
                route.native_table(reference, patients, ['A', 'B'], ['sample'])


class Observer(unittest.TestCase):
    def test_observer_preserves_script_and_records_native_return_values(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            script = root / 'native.py'
            script.write_text('''import numpy as np
from types import SimpleNamespace
args = SimpleNamespace(random_restarts=2)
samples = ['sample']
tissues = ['A', 'B']
def em(seed):
    i = 3
    convergence_criteria = .0001
    alpha_diff, gamma_diff = .00001, .00002
    return np.array([[.7, .3]]), np.array([[.5], [.5]]), float(seed)
for seed in (7, 11):
    em(seed)
''')
            original = sha(script)
            before = sys.getprofile()
            result = observe(script, root / 'diagnostics.json', [])
            self.assertEqual(sys.getprofile(), before)
            self.assertEqual(sha(script), original)
            self.assertEqual([r['seed'] for r in result['restarts']], [7, 11])
            self.assertEqual(result['best_restart'], 1)
            self.assertTrue(result['all_starts_converged'])
            self.assertEqual(result['restarts'][0]['alpha'], [[.7, .3]])

    def test_observer_retains_native_exception_and_restores_profile(self):
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / 'native.py'
            script.write_text("def em():\n    raise RuntimeError('native failure')\nem()\n")
            before = sys.getprofile()
            with self.assertRaisesRegex(RuntimeError, 'native failure'):
                observe(script, Path(directory) / 'diagnostics.json', [])
            self.assertEqual(sys.getprofile(), before)


class Execution(unittest.TestCase):
    def setUp(self):
        guard = patch.object(route.output_guard, 'require_clean_source', return_value='a' * 40)
        guard.start(); self.addCleanup(guard.stop)
        for name in ('resume_nextflow_session', 'save_nextflow_session'):
            patcher = patch.object(route, name, return_value='01234567-1234-1234-1234-123456789abc')
            patcher.start(); self.addCleanup(patcher.stop)


    def test_tool_specific_version_commands(self):
        with patch.object(route.subprocess, 'check_output', return_value='pinned version') as call:
            self.assertEqual(route.runtime_versions(), dict(nextflow='pinned version', apptainer='pinned version'))
        self.assertEqual([item.args[0] for item in call.call_args_list],
                         [['nextflow', '-version'], ['apptainer', '--version']])

    def test_failed_preflight_records_start_and_preserves_it_on_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / '2026-10-26_12-00-celfie-test-v1'
            argv = ['celfie_patients.py', str(output), '--runtime', str(Path(directory) / 'missing-runtime'),
                    '--reference-counts', 'reference.tsv', '--patient-counts', 'patients.tsv',
                    '--reference-validation', 'gate.json', '--samples', 'sample']

            def create(path):
                route.output_guard.reserve_path(path, {'project_revision': 'a' * 40, 'code_sha256': sha(route.WORKFLOW/'common.py')})
                route.output_guard.update(path,status='running')
                write_json(path / 'execution.json', dict(started_utc='2026-10-26T12:00:00+00:00', status='running'))
                return path

            with patch.dict(os.environ, SLURM_JOB_ID='mock'), \
                 patch.object(route, 'new_run', side_effect=create) as new_run, \
                 patch('shutil.which', return_value='/mock'), \
                 patch.object(route, 'verify_runtime', side_effect=RuntimeError('preflight runtime failure')), \
                 patch.object(route, 'invoke') as invoke:
                for resume in (False, True):
                    with patch.object(sys, 'argv', argv + (['--resume'] if resume else [])):
                        with self.assertRaisesRegex(RuntimeError, 'preflight runtime failure'):
                            route.main()
                    state = read_json(output / 'execution.json')
                    self.assertEqual(state['started_utc'], '2026-10-26T12:00:00+00:00')
                    self.assertEqual(state['status'], 'failed')
                    self.assertEqual(state['stage'], 'preflight')
                self.assertEqual(new_run.call_count, 1)
                invoke.assert_not_called()
                records = [json.loads(line) for line in (output / 'preflight-attempts.jsonl').read_text().splitlines()]
                self.assertEqual([(r['attempt'], r['status']) for r in records],
                                 [(1, 'running'), (1, 'failed'), (2, 'running'), (2, 'failed')])

    def test_native_trace_container_and_published_output_are_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            runtime = root / 'runtime'
            task = root / 'task'
            task.mkdir()
            (root / 'results').mkdir()
            fractions = '\tA\tB\nsample\t0.7\t0.3\n'
            (task / '1_tissue_proportions.txt').write_text(fractions)
            (root / 'results/1_tissue_proportions.txt').write_text(fractions)
            for name in ('.command.sh', '.command.out', '.command.err'):
                (task / name).write_text('native evidence\n')
            (task / '.command.run').write_text(f'apptainer exec {runtime}/images/python.sif native\n')

            def native(command, *_args, **_kwargs):
                trace = Path(command[command.index('--trace_path') + 1])
                trace.write_text(f'name\tstatus\tworkdir\nRUN_DECONV\tCACHED\t{task}\n')

            route.output_guard._write(root, dict(schema=1, run_id=root.name, output=str(root.resolve()),
                identity={'project_revision': 'a' * 40, 'code_sha256': sha(route.WORKFLOW/'common.py')},status='running',
                started_utc='2026-11-02T13:00:00Z',session_id=None,attempts=[]))
            with patch.object(route, 'run', side_effect=native) as run:
                route.invoke(root, runtime, dict(diagnostics=False), 'fit', 1,
                             resume=True, require_cache=True)
                self.assertEqual(run.call_args.args[0][1], '-C')
                self.assertTrue((root / 'native-task-evidence/fit.command.run').exists())
                (root / 'results/1_tissue_proportions.txt').write_text('stale fractions\n')
                with self.assertRaisesRegex(ValueError, 'Published fractions differ'):
                    route.invoke(root, runtime, dict(diagnostics=False), 'recovery', 1,
                                 resume=True, require_cache=True)

    def test_mocked_fit_cached_recovery_and_config_drift_boundary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            runtime = root / 'runtime'
            runtime.mkdir()
            write_json(runtime / 'runtime.lock.json', {})
            counts(root / 'reference.tsv', ['A', 'B'])
            counts(root / 'patients.tsv', ['sample'])
            write_json(root / 'gate.json', dict(status='passed', atlas_sha256=ATLAS_SHA,
                       reference_counts_sha256=sha(root / 'reference.tsv')))
            output = root / '2026-10-26_12-00-celfie-test-v1'
            argv = ['celfie_patients.py', str(output), '--runtime', str(runtime),
                    '--reference-counts', str(root / 'reference.tsv'),
                    '--patient-counts', str(root / 'patients.tsv'),
                    '--reference-validation', str(root / 'gate.json'),
                    '--samples', 'sample', '--verify-cache']

            def create(path):
                route.output_guard.reserve_path(path, {'project_revision': 'a' * 40, 'code_sha256': sha(route.WORKFLOW/'common.py')})
                route.output_guard.update(path,status='running')
                write_json(path / 'execution.json', dict(started_utc='2026-10-26T12:00:00+00:00', status='running'))
                return path

            def fit(path, _runtime, _params, _phase, _attempt, **_kwargs):
                (path / 'results').mkdir(exist_ok=True)
                (path / 'results/1_tissue_proportions.txt').write_text('\tA\tB\nsample\t0.7\t0.3\n')

            with patch.dict(os.environ, SLURM_JOB_ID='mock'), \
                 patch.object(route, 'verify_runtime', return_value=runtime), \
                 patch.object(route, 'atlas_rows', return_value=([
                     dict(chr='chr1', start='10', end='20'), dict(chr='chr2', start='30', end='40')], ['A', 'B'])), \
                 patch.object(route, 'new_run', side_effect=create), \
                 patch('shutil.which', return_value='/mock'), \
                 patch.object(route.subprocess, 'check_output', return_value='pinned version'), \
                 patch.object(route, 'invoke', side_effect=fit) as invoke:
                with patch.object(sys, 'argv', argv):
                    route.main()
                self.assertEqual(invoke.call_count, 2)
                self.assertTrue(read_json(output / 'COMPLETE.json')['cached_recovery_verified'])
                original_start = read_json(output / 'execution.json')['started_utc']
                (output / 'celfie.config').write_text("process.executor = 'local'\n")
                invoke.reset_mock()
                with patch.object(sys, 'argv', argv + ['--resume']), self.assertRaisesRegex(ValueError, 'incomplete|immutable'):
                    route.main()
                invoke.assert_not_called()
                self.assertEqual(read_json(output / 'execution.json')['started_utc'], original_start)
                attempts = [json.loads(line) for line in (output / 'attempts.jsonl').read_text().splitlines()]
                self.assertEqual([(a['attempt'], a['status']) for a in attempts],
                                 [(1, 'running'), (1, 'complete')])


if __name__ == '__main__':
    unittest.main()
