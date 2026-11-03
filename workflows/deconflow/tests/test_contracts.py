"""Scientific boundary tests that must pass before native execution."""
import importlib.util
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.reference import equal_sample_beta, membership
from native import checked
from native import invoke
from tools.compare import differences, compare
from analysis import selected_samples
from common import new_run, table, write_json, sha, output_guard, WORKFLOW
from fragments import prepare_fragments, compare_fragments


def _own_fixture(root):
    if not (root / output_guard.RECORD).exists():
        output_guard._write(root, dict(schema=1, run_id=root.name, output=str(root.resolve()),
            identity={'project_revision': 'a' * 40, 'code_sha256': sha(WORKFLOW/'common.py')}, status='running',
            started_utc='2026-11-02T13:00:00Z', session_id=None, attempts=[]))


class Contracts(unittest.TestCase):
    def setUp(self):
        guard = patch.object(output_guard, 'require_clean_source', return_value='a' * 40)
        guard.start(); self.addCleanup(guard.stop)
        # Low-level invocation tests isolate subprocess execution; session parsing
        # and cache/history identity are exercised in test_nextflow_sessions.py.
        for target, value in (('native.resume_nextflow_session', '01234567-1234-1234-1234-123456789abc'),
                              ('native.save_nextflow_session', '01234567-1234-1234-1234-123456789abc')):
            patcher = patch(target, return_value=value)
            patcher.start(); self.addCleanup(patcher.stop)


    def test_fragment_triplet_comparison_selects_from_full_saved_cohort(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            records = [dict(sample=f'{p}-{c}', patient=p, compartment=c) for p in
                       ['p1', 'p2', 'p3', 'p4'] for c in ['PB', 'BM', 'cfDNA']]
            names = [r['sample'] for r in records]
            table(root / 'samples.tsv', records, ['sample', 'patient', 'compartment'])
            (root / 'inputs').mkdir()
            baseline = root / 'saved'
            (baseline / 'results').mkdir(parents=True)
            saved = pd.DataFrame([[.6, .4]] * 12, index=names, columns=['A', 'B'])
            saved.T.to_csv(baseline / 'results/lineage_fractions.csv')
            saved.T.to_csv(baseline / 'epidish_rpc_fractions.tsv', sep='\t')
            write_json(root / 'baselines.json', [dict(cohort='test', uxm=str(baseline))])
            write_json(root / 'inputs/triplet.jobs.json', [dict(id='test_ge1', folder=str(baseline),
                       samples=names[:3], cohort='test', threshold=1)])
            for method in ('uxm', 'epidish'):
                dest = root / 'triplet' / method
                if method == 'epidish':
                    dest /= 'test_ge1'
                dest /= 'deconvolution_results/reference_based/' + method
                dest.mkdir(parents=True)
                saved.loc[names[:3]].to_csv(dest / ('output.csv' if method == 'uxm' else 'epidish_res_RPC.csv'))
            result = compare_fragments(root, 'triplet', ['A', 'B'])
            self.assertEqual(result['uxm']['estimates'], 6)
            self.assertEqual(result['epidish']['maximum_absolute_difference_pp'], 0)

    def test_fragment_adapter_preserves_cohort_mask_and_values(self):
        from common import read_json
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            saved = root / 'saved'
            saved.mkdir()
            ids = ['chr1:10-20|U', 'chr2:30-40|M']
            records = [dict(sample=f'{p}-{c}', patient=p, compartment=c) for p in
                       ['p1', 'p2', 'p3', 'p4'] for c in ['PB', 'BM', 'cfDNA']]
            names = [r['sample'] for r in records]
            ref = pd.DataFrame([[.9, .1], [.2, .8]], index=ids, columns=['A', 'B'])
            ref.to_csv(saved / 'reference_matrix.tsv', sep='\t')
            pd.DataFrame([[.3] * 12, [.7] * 12], index=ids, columns=names).to_csv(
                saved / 'sample_regional_values.tsv.gz', sep='\t')
            pd.DataFrame(dict(chr=['chr1', 'chr2'], start=[10, 30], end=[20, 40],
                              name=['chr1:10-20', 'chr2:30-40'])).to_csv(saved / 'marker_mask.tsv', sep='\t', index=False)
            artifact = root / 'existing.pat'
            artifact.write_text('fixture')
            inventory = dict(pats={s: {k: {'path': str(artifact)} for k in ('pat', 'csi')} for s in names},
                             jobs=[dict(id='test_ge1', folder=str(saved), samples=names)])
            prepare_fragments(root, records, inventory)
            triplet = read_json(root / 'inputs/triplet.jobs.json')[0]
            actual = pd.read_csv(triplet['matrix'], index_col=0)
            reference = pd.read_csv(triplet['reference'], index_col=0)
            self.assertEqual(list(actual.columns), ['chr', 'start', 'end'] + names[:3])
            self.assertEqual(list(actual.index), list(reference.index))
            np.testing.assert_allclose(reference.values, ref.values)
            np.testing.assert_allclose(actual[names[:3]].values, [[.3] * 3, [.7] * 3])
            self.assertEqual(len(actual), 2)  # Selecting samples does not recalculate the cohort mask.
            self.assertTrue((root / 'inputs' / (names[0] + '.pat.gz')).is_symlink())

    def test_patient_entrypoint_does_not_import_migration_tools(self):
        code = "import analysis, sys; assert not any(n.startswith(('tools.', 'matplotlib')) for n in sys.modules)"
        subprocess.run([sys.executable, '-c', code], check=True, cwd=Path(__file__).resolve().parents[1])

    def test_native_invocations_keep_settings_and_resume(self):
        from common import read_json, WORK
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for method, sheet_key, atlas_key in [('uxm', 'test_bams', 'uxm_atlas'), ('epidish', 'test_set', 'ref_matrix')]:
                result = root / 'triplet' / method
                (result / 'pipeline_info').mkdir(parents=True)
                table(result / 'pipeline_info/execution_trace_triplet.txt',
                      [dict(status='COMPLETED')], ['status'])
                with patch('native.run') as run:
                    _own_fixture(root)
                    invoke(root / 'runtime', root, method, root / 'samples.csv', root / 'atlas.tsv', 'triplet')
                    params = read_json(root / 'native' / method / 'triplet.params.json')
                    self.assertEqual(params[sheet_key], str(root / 'samples.csv'))
                    self.assertEqual(params[atlas_key], str(root / 'atlas.tsv'))
                    self.assertEqual({k: params[k] for k in ('genome', 'DMRselection', 'benchmark', 'bulk_min_counts', 'bulk_min_cpgs')},
                                     dict(genome='hg38', DMRselection=None, benchmark=False, bulk_min_counts=0, bulk_min_cpgs=0))
                    _own_fixture(root)
                    invoke(root / 'runtime', root, method, root / 'samples.csv', root / 'atlas.tsv', 'triplet', resume=True)
                    command = run.call_args.args[0]
                    self.assertIn(WORK, command)
                    self.assertEqual(command[-2:], ['-resume', '01234567-1234-1234-1234-123456789abc'])
                    self.assertEqual(run.call_args.kwargs['env']['DECONFLOW_RUNTIME'], str(root / 'runtime'))
                    with self.assertRaises(ValueError):
                        invoke(root / 'runtime', root, method, root / 'samples.csv', root / 'changed.tsv', 'triplet', resume=True)

    def test_fragment_invocation_uses_native_modules_and_separate_sessions(self):
        from common import WORKFLOW, read_json, COMMIT
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for tag in ('cohort_full', 'cohort_hematopoietic'):
                result = root / tag / 'uxm'
                (result / 'pipeline_info').mkdir(parents=True)
                table(result / 'pipeline_info' / f'execution_trace_{tag}.txt',
                      [dict(status='COMPLETED')], ['status'])
                with patch('native.run') as run:
                    _own_fixture(root)
                    invoke(root / 'runtime', root, 'uxm', root / 'samples.csv', root / 'atlas.tsv',
                           tag, fragments=True)
                    command = run.call_args.args[0]
                    self.assertIn(WORKFLOW / 'fragments.nf', command)
                    self.assertNotIn('-profile', command)  # Enabled explicitly by Flamingo config.
                    self.assertEqual(run.call_args.kwargs['cwd'], root / 'native/uxm' / tag)
                    params = read_json(root / 'native/uxm' / tag / (tag + '.params.json'))
                    self.assertEqual(params['deconflow_source'], str(root / 'runtime' / ('DecoNFlow-' + COMMIT)))
                    self.assertNotIn('maxit', params)

    def test_local_pat_entrypoint_rejects_epidish_before_launch(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            with patch('native.run') as run:
                with self.assertRaisesRegex(ValueError, 'EpiDISH is retired'):
                    _own_fixture(root)
                    invoke(root / 'runtime', root, 'epidish', root / 'samples.csv',
                           root / 'atlas.csv', 'cohort', fragments=True)
                run.assert_not_called()
            self.assertFalse((root / 'native').exists())

    def test_pat_plan_runs_only_native_uxm_and_records_method(self):
        import os
        from types import SimpleNamespace
        from analysis import run_pat_plan
        from common import read_json
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            runtime = root / 'runtime'
            runtime.mkdir()
            write_json(runtime / 'runtime.lock.json', {})
            output = root / '2026-10-26_10-00-pat-fixture-v1'
            inventory = dict(records=[dict(sample='new_sample')],
                             settings=dict(profiles=['full', 'hematopoietic']))
            args = SimpleNamespace(output=output, phase='cohort', plan=root / 'plan.json', resume=False)

            def create_run(path):
                output_guard.reserve_path(path, {'project_revision': 'a' * 40, 'code_sha256': sha(WORKFLOW/'common.py')})
                output_guard.update(path,status='running')
                write_json(path / 'execution.json', dict(started_utc='2026-10-26T10:00:00+00:00'))
                return path

            def prepare(path, plan):
                (path / 'inputs').mkdir()
                (path / 'inputs/cohort.pats.csv').write_text('sample,pat,csi\n')

            with patch.dict(os.environ, SLURM_JOB_ID='fixture'), \
                    patch('analysis.read_plan', return_value=inventory), \
                    patch('analysis.new_run', side_effect=create_run), \
                    patch('analysis.prepare_pat_inputs', side_effect=prepare), \
                    patch('analysis.profile_atlas', return_value=(None, ['A', 'B'])), \
                    patch('analysis.invoke') as run, patch('analysis.native_fractions') as fractions:
                run_pat_plan(args, runtime)
            self.assertEqual([call.args[2] for call in run.call_args_list], ['uxm', 'uxm'])
            self.assertEqual([call.args[5] for call in run.call_args_list],
                             ['cohort_full', 'cohort_hematopoietic'])
            self.assertEqual(fractions.call_count, 2)
            self.assertEqual(read_json(output / 'execution.json')['methods'], ['uxm'])
            self.assertTrue((output / 'COHORT_COMPLETE').is_file())
            self.assertFalse((output / 'inputs/prepared').exists())

    def test_equal_sample_weighting_and_missing(self):
        # An observed zero contributes; absent coverage does not. Pooled counts
        # would give 0.9 here, while the approved equal-sample mean is 0.5.
        counts = np.array([[[0, 1], [0, 0]], [[9, 9], [3, 4]], [[2, 4], [0, 4]]])
        beta, n = equal_sample_beta(counts, ['A', 'A', 'B'], ['A', 'B'])
        np.testing.assert_allclose(beta, [[.5, .5], [.75, 0]])
        np.testing.assert_array_equal(n, [[2, 1], [1, 1]])

    def test_missing_reference_cell_or_marker_stops(self):
        with self.assertRaises(ValueError):
            equal_sample_beta(np.array([[[0, 0]]]), ['A'], ['A'])
        with self.assertRaises(ValueError):
            equal_sample_beta(np.array([[[1, 1]]]), ['A'], ['A', 'B'])

    def test_invalid_counts_stop(self):
        with self.assertRaises(ValueError):
            equal_sample_beta(np.array([[[2, 1]]]), ['A'], ['A'])

    def test_maintained_wgbstools_cpg_boundaries(self):
        from common import PROJECT
        tool_source = Path.home() / '.conda/envs/nextflow/share/wgbs_tools/src/python'
        sys.path.insert(0, str(tool_source))
        from beta_to_blocks import collapse_process
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'counts.lbeta'
            np.array([[1, 2], [0, 0], [3, 4], [2, 2]], dtype=np.uint16).tofile(p)
            blocks = pd.DataFrame([['chr1', 10, 20, 1, 3], ['chr1', 20, 30, 3, 5]],
                                  columns=['chr', 'start', 'end', 'startCpG', 'endCpG'])
            np.testing.assert_array_equal(collapse_process(str(p), blocks, False), [[1, 2], [5, 6]])

    def test_full_atlas_membership_gate(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'mapping.tsv'
            fields = ['gsm_id', 'status', 'atlas_cell_type', 'pat', 'evidence_url']
            table(path, [dict(gsm_id='GSM1', status='verified', atlas_cell_type='A', pat='/a', evidence_url='source')], fields)
            with self.assertRaises(ValueError):
                membership(path, ['A', 'Megakaryocytes'])
            table(path, [dict(gsm_id='GSM1', status='unresolved', atlas_cell_type='A', pat='/a', evidence_url='source')], fields)
            with self.assertRaises(ValueError):
                membership(path, ['A'])

    def test_labels_and_difference_sign(self):
        a = pd.DataFrame([[.2, .8]], index=['sample'], columns=['B', 'A'])
        a = checked(a, ['sample'], ['A', 'B'])
        b = pd.DataFrame([[.7, .3]], index=['sample'], columns=['A', 'B'])
        result = differences(a, b, 'uxm', 'saved')
        self.assertAlmostEqual(result[0]['difference_pp'], 10)
        self.assertAlmostEqual(result[1]['difference_pp'], -10)
        with self.assertRaises(ValueError):
            checked(a, ['other'], ['A', 'B'])
        with self.assertRaises(ValueError):
            checked(a * 2, ['sample'], ['A', 'B'])

    def test_patient_triplet_selection(self):
        data = [dict(sample=f'{p}-{c}', patient=p, compartment=c) for p in ['p2', 'p1', 'p3', 'p4'] for c in ['PB', 'BM', 'cfDNA']]
        self.assertEqual([r['sample'] for r in selected_samples(data, 'triplet')], ['p2-PB', 'p2-BM', 'p2-cfDNA'])
        data[-1]['compartment'] = 'BM'
        with self.assertRaises(ValueError):
            selected_samples(data, 'cohort')

    def test_invalid_datetime_and_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                new_run(Path(d) / '2026-02-30_14-00-example-v1')
            with self.assertRaises(ValueError):
                new_run(Path(d) / '2026-02-28_25-00-example-v1')
            p = Path(d) / '2026-02-28_14-00-example-v1'
            new_run(p)
            with self.assertRaises(FileExistsError):
                new_run(p)

    def test_saved_comparison_export_and_checksum_gate(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            analysis = root / 'analysis'
            analysis.mkdir()
            (analysis / 'COHORT_COMPLETE').write_text('synthetic')
            write_json(analysis / 'input_fingerprint.json', {'test': True})
            data = [dict(sample=f'{p}-{c}', patient=p, compartment=c) for p in ['p1', 'p2', 'p3', 'p4'] for c in ['PB', 'BM', 'cfDNA']]
            samples = [r['sample'] for r in data]
            cells = ['A', 'B']
            table(analysis / 'samples.tsv', data, ['sample', 'patient', 'compartment'])
            frame = pd.DataFrame([[.7, .3]] * 12, index=samples, columns=cells)
            for method in ('uxm', 'epidish'):
                folder = analysis / 'cohort' / method
                dest = folder / 'deconvolution_results/reference_based' / method
                dest.mkdir(parents=True)
                if method == 'uxm':
                    frame.to_csv(dest / 'output.csv')
                else:
                    frame.to_csv(dest / 'epidish_res_RPC.csv')
                (folder / 'pipeline_info').mkdir()
                table(folder / 'pipeline_info/execution_trace_cohort.txt', [], ['name', 'workdir'])
            baseline = root / 'baseline'
            (baseline / 'results').mkdir(parents=True)
            frame.T.to_csv(baseline / 'results/lineage_fractions.csv')
            (baseline / 'provenance').mkdir()
            write_json(baseline / 'provenance/input_manifest.json', {'epidish_maxit': 75,
                       'full_uxm_run': {'native_result_sha256': sha(baseline / 'results/lineage_fractions.csv')}})
            write_json(analysis / 'baselines.json', [{'cohort': 'synthetic', 'uxm': str(baseline), 'epidish': str(baseline)}])
            for threshold in (1, 3):
                folder = baseline / f'comparisons/full/coverage-ge{threshold}'
                folder.mkdir(parents=True)
                frame.T.to_csv(folder / 'epidish_rpc_fractions.tsv', sep='\t')
                (folder / 'marker_mask.tsv').write_text('marker\nexample\n')
                (folder / 'epidish_warnings.tsv').write_text('sample\twarning\n')
                write_json(folder / 'manifest.json', {'profile': 'full', 'coverage_threshold': threshold,
                           'samples': samples, 'marker_count': 1,
                           'files': [{'name': name, 'sha256': sha(folder / name)} for name in
                                     ['epidish_rpc_fractions.tsv', 'epidish_warnings.tsv', 'marker_mask.tsv']]})
            result = new_run(root / '2026-10-05_00-01-synthetic-comparison-v1')
            with patch('tools.compare.atlas_rows', return_value=([], cells)):
                compare(result, analysis, 'cohort')
                out = pd.read_csv(result / 'all_differences.tsv', sep='\t')
                self.assertEqual(len(out), 12 * 2 * 3)
                self.assertTrue((out.difference_pp == 0).all())
                self.assertTrue((result / 'scatterplots.pdf').is_file())
                (baseline / 'comparisons/full/coverage-ge1/marker_mask.tsv').write_text('changed')
                other = new_run(root / '2026-10-05_00-02-synthetic-comparison-v1')
                with self.assertRaises(ValueError):
                    compare(other, analysis, 'cohort')


if __name__ == '__main__':
    unittest.main()
