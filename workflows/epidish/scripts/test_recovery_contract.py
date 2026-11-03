#!/usr/bin/env python3
"""Temporary EpiDISH contract fixtures; never run a scientific workflow."""
import copy
import csv
import gzip
import json
import os
import shlex
import subprocess
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import recovery_contract as contract

SCRIPTS = Path(__file__).resolve().parent


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + '\n')


def put(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == '.gz':
        with gzip.open(path, 'wt') as handle:
            handle.write(content)
    else:
        path.write_text(content)


def files_record(directory, names):
    return [{'name': name, 'bytes': (directory / name).stat().st_size,
             'sha256': contract.sha256_file(directory / name)} for name in names]


class ContractFixtures(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='epidish-contract-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.work = self.root / 'work'
        self.baseline = self.root / 'relocated-baseline'
        self.inputs = {
            'samples': ['s1', 's2'], 'profiles': ['full', 'hematopoietic'],
            'coverage_thresholds': [1, 3, 5], 'regional_measure': 'directional U or M / total',
            'epidish_maxit': 50, 'atlas': {'sha256': 'a' * 64},
            'uxm_checkout': {'commit': 'b' * 40},
            'pat_artifacts': [{'sample': sample, 'kind': kind, 'bytes': 7,
                               'sha256': str(index + 1) * 64,
                               'path': '/historical/path', 'modified_utc': 'historical'}
                              for index, (sample, kind) in enumerate(
                                  (s, k) for s in ('s1', 's2') for k in ('pat', 'csi', 'lbeta'))],
        }
        fractions = 'CellType\ts1\ts2\nBlood-B\t0.7\t0.2\nBlood-T\t0.3\t0.8\n'
        put(self.root / 'native.csv', 'CellType,s1,s2\nBlood-B,0.7,0.2\nBlood-T,0.3,0.8\n')
        self.jobs = []
        self.comparisons = []
        for profile in ('full', 'hematopoietic'):
            for threshold in (1, 3, 5):
                label = f'coverage-ge{threshold}'
                relative = f'comparisons/{profile}/{label}'
                contents = {
                    'marker_mask.tsv': 'region_id\tminimum_cohort_coverage\tcoverage_s1\tcoverage_s2\nr1\t10\t10\t10\nr2\t10\t10\t10\nr3\t10\t10\t10\n',
                    'coverage_summary.tsv': 'sample\tmarker_count\ns1\t3\ns2\t3\n',
                    'sample_regional_values.tsv.gz': 'region_id\ts1\ts2\nr1\t0.34\t0.74\nr2\t0.5\t0.5\nr3\t0.66\t0.26\n',
                    'reference_matrix.tsv': 'region_id\tBlood-B\tBlood-T\nr1\t0.1\t0.9\nr2\t0.5\t0.5\nr3\t0.9\t0.1\n',
                    'epidish_rpc_fractions.tsv': fractions,
                    'epidish_warnings.tsv': 'sample\twarning\tprofile\tcoverage_threshold\n',
                    'uxm_native_fractions.tsv': fractions,
                    'uxm_matched_atlas.l4.tsv': 'region\tBlood-B\tBlood-T\nr1\t0.1\t0.9\n',
                }
                for name, content in contents.items():
                    put(self.work / relative / name, content)
                    put(self.baseline / relative / name, content)
                put(self.baseline / relative / 'uxm_weighted_matched_fractions.tsv', fractions)
                comparison = {'profile': profile, 'coverage_threshold': threshold, 'marker_count': 3,
                              'sample_count': 2, 'samples': ['s1', 's2'],
                              'files': files_record(self.baseline / relative,
                                                    [*contract.PREPARED_FILES, 'uxm_weighted_matched_fractions.tsv'])}
                self.comparisons.append(comparison)
                write_json(self.baseline / relative / 'manifest.json', comparison)
                self.jobs.append([profile, threshold, label, 3, relative,
                                  f'{relative}/uxm_matched_atlas.l4.tsv', str(self.root / 'native.csv')])
        put(self.work / 'provenance/comparison_jobs.tsv',
            '\t'.join(['profile', 'coverage_threshold', 'threshold_label', 'marker_count',
                       'comparison_dir', 'matched_atlas', 'native_source']) + '\n' +
            ''.join('\t'.join(map(str, job)) + '\n' for job in self.jobs))
        put(self.work / 'provenance/pat_manifest.tsv', 'sample\tpat\tcsi\tlbeta\ns1\tp1\tc1\tl1\ns2\tp2\tc2\tl2\n')
        put(self.work / 'provenance/sample_metadata.tsv', 'sample\ts1\ns1\tvalue\n')
        write_json(self.work / 'provenance/recovery_contract.json', {'frozen': True})
        write_json(self.work / 'provenance/execution_attempts.json', {
            'started_utc': '2026-10-28T08:17:21Z',
            'attempts': [{'started_utc': '2026-10-28T08:17:21Z', 'attempt_id': '1', 'status': 'failed'}]})
        write_json(self.work / 'provenance/input_manifest.json', self.inputs)
        write_json(self.baseline / 'provenance/input_manifest.json', self.inputs)
        put(self.baseline / 'provenance/commands.log', 'historical command\n')
        put(self.baseline / 'comparisons/summary/marker_counts.tsv', 'profile\tcount\nfull\t3\n')
        put(self.baseline / 'EPIDISHR_COMPLETED', 'completed\n')
        self.refresh_manifest()
        write_json(self.work / 'provenance/recovery_contract.json', {
            'frozen': True, 'baseline': {'directory': str(self.baseline),
                'run_manifest': contract.record(self.baseline / 'provenance/run_manifest.json')}})
        for job in self.jobs:
            contract.write_receipt(self.work, f'prepared/{job[0]}/{job[2]}',
                                   [f'{job[4]}/{name}' for name in contract.PREPARED_FILES])
        contract.write_receipt(self.work, 'preparation', [
            'provenance/comparison_jobs.tsv', 'provenance/pat_manifest.tsv',
            'provenance/sample_metadata.tsv', 'provenance/input_manifest.json'])

    def refresh_manifest(self):
        manifest = {'completed_utc': '2026-09-11T12:34:56Z', 'samples': ['s1', 's2'],
                    'comparison_count': 6, 'comparisons': self.comparisons,
                    'input_manifest_sha256': contract.sha256_file(self.baseline / 'provenance/input_manifest.json'),
                    'commands_log_sha256': contract.sha256_file(self.baseline / 'provenance/commands.log'),
                    'summary_files': files_record(self.baseline / 'comparisons/summary', ['marker_counts.tsv'])}
        write_json(self.baseline / 'provenance/run_manifest.json', manifest)

    def run_reuse(self, success=True):
        result = subprocess.run(['Rscript', str(SCRIPTS / 'reuse_matched_uxm.R'),
                                 str(self.work), str(self.baseline)], capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(len(list(self.work.glob('comparisons/*/*/.uxm_weighted_matched.raw.csv'))), 6)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(list(self.work.glob('comparisons/*/*/.uxm_weighted_matched.raw.csv')), [], result.stderr)
        return result

    def test_compatible_relocated_baseline_and_rpc_only_sensitivity(self):
        changed = copy.deepcopy(self.inputs)
        changed['epidish_maxit'] = 200
        for record in changed['pat_artifacts']:
            record['path'] = '/relocated/current/path'
            record['modified_utc'] = 'new timestamp'
        write_json(self.work / 'provenance/input_manifest.json', changed)
        contract.receipt_path(self.work, 'preparation').unlink()
        contract.write_receipt(self.work, 'preparation', [
            'provenance/comparison_jobs.tsv', 'provenance/pat_manifest.tsv',
            'provenance/sample_metadata.tsv', 'provenance/input_manifest.json'])
        self.run_reuse()
        for job in self.jobs:
            self.assertTrue(contract.check_receipt(self.work, f'matched/{job[0]}/{job[2]}',
                                                   [f'{job[4]}/.uxm_weighted_matched.raw.csv']))
        self.run_reuse()  # recovery revalidates rather than rewriting sealed reuse provenance

    def test_changed_scientific_or_tool_identity_fails_before_import(self):
        for field, value in [('atlas', {'sha256': 'f' * 64}), ('uxm_checkout', {'commit': 'f' * 40}),
                             ('regional_measure', 'different'), ('samples', ['s2', 's1']),
                             ('profiles', ['hematopoietic', 'full']), ('coverage_thresholds', [1, 5, 10])]:
            with self.subTest(field=field):
                changed = copy.deepcopy(self.inputs)
                changed[field] = value
                write_json(self.baseline / 'provenance/input_manifest.json', changed)
                self.refresh_manifest()
                self.run_reuse(False)

    def test_missing_hash_or_size_evidence_fails(self):
        changed = copy.deepcopy(self.inputs)
        del changed['pat_artifacts'][-1]['bytes']
        write_json(self.baseline / 'provenance/input_manifest.json', changed)
        self.refresh_manifest()
        self.run_reuse(False)

    def test_changed_payload_identity_fails(self):
        changed = copy.deepcopy(self.inputs)
        changed['pat_artifacts'][-1]['sha256'] = 'f' * 64
        write_json(self.baseline / 'provenance/input_manifest.json', changed)
        self.refresh_manifest()
        self.run_reuse(False)

    def test_late_mask_difference_fails_before_any_import(self):
        job = self.jobs[-1]
        directory = self.work / job[4]
        put(directory / 'marker_mask.tsv', (directory / 'marker_mask.tsv').read_text().replace('r3', 'changed'))
        stage = f'prepared/{job[0]}/{job[2]}'
        contract.receipt_path(self.work, stage).unlink()
        contract.write_receipt(self.work, stage, [f'{job[4]}/{name}' for name in contract.PREPARED_FILES])
        self.run_reuse(False)

    def test_reference_or_sample_values_difference_fails(self):
        for name in ('reference_matrix.tsv', 'sample_regional_values.tsv.gz'):
            with self.subTest(name=name):
                path = self.work / self.jobs[-1][4] / name
                original = path.read_bytes()
                if path.suffix == '.gz':
                    with gzip.open(path, 'rt') as stream:
                        text = stream.read()
                else:
                    text = path.read_text()
                put(path, text.replace('0.', '0.0', 1))
                job = self.jobs[-1]
                stage = f'prepared/{job[0]}/{job[2]}'
                contract.receipt_path(self.work, stage).unlink()
                contract.write_receipt(self.work, stage, [f'{job[4]}/{n}' for n in contract.PREPARED_FILES])
                self.run_reuse(False)
                path.write_bytes(original)
                contract.receipt_path(self.work, stage).unlink()
                contract.write_receipt(self.work, stage, [f'{job[4]}/{n}' for n in contract.PREPARED_FILES])

    def test_tampered_fractions_fails_before_any_import(self):
        path = self.baseline / self.jobs[-1][4] / 'uxm_weighted_matched_fractions.tsv'
        path.write_text(path.read_text().replace('0.7', '0.6'))
        self.run_reuse(False)

    def test_partial_matched_output_without_receipt_is_refused(self):
        path = self.work / self.jobs[0][4] / '.uxm_weighted_matched.raw.csv'
        put(path, 'partial\n')
        with self.assertRaisesRegex(ValueError, 'Unreceipted'):
            contract.check_recovery_outputs(self.work)
        self.run_reuse(False) if not path.exists() else self.assertEqual(path.read_text(), 'partial\n')

    def test_missing_legacy_contract_fails_without_output_mutation(self):
        legacy = self.root / 'legacy'
        legacy.mkdir()
        put(legacy / 'unchanged', 'original\n')
        before = {p.relative_to(legacy): p.read_bytes() for p in legacy.rglob('*') if p.is_file()}
        with self.assertRaisesRegex(ValueError, 'Missing recovery contract'):
            contract.reserve(Namespace(work_dir=legacy, recover=True))
        self.assertEqual(before, {p.relative_to(legacy): p.read_bytes() for p in legacy.rglob('*') if p.is_file()})

    def test_valid_recovery_preserves_original_contract_and_timing(self):
        original_contract = (self.work / 'provenance/recovery_contract.json').read_bytes()
        with patch.object(contract, 'build_contract', return_value=json.loads((self.work / 'provenance/recovery_contract.json').read_text())):
            contract.reserve(Namespace(work_dir=self.work, recover=True,
                                       started_utc='2026-10-28T10:25:31Z', attempt_id='2'))
        self.assertEqual(original_contract, (self.work / 'provenance/recovery_contract.json').read_bytes())
        timing = json.loads((self.work / 'provenance/execution_attempts.json').read_text())
        self.assertEqual(timing['started_utc'], '2026-10-28T08:17:21Z')
        self.assertEqual(len(timing['attempts']), 2)
        contract.finish(Namespace(work_dir=self.work, attempt_id='2', status='failed'))
        timing = json.loads((self.work / 'provenance/execution_attempts.json').read_text())
        self.assertNotIn('completed_utc', timing['attempts'][-1])
        self.assertIn('ended_utc', timing['attempts'][-1])

    def test_contract_change_and_unreceipted_preparation_fail_before_attempt_append(self):
        timing = (self.work / 'provenance/execution_attempts.json').read_bytes()
        with patch.object(contract, 'build_contract', return_value={'frozen': False}):
            with self.assertRaisesRegex(ValueError, 'contract differs'):
                contract.reserve(Namespace(work_dir=self.work, recover=True))
        self.assertEqual(timing, (self.work / 'provenance/execution_attempts.json').read_bytes())
        contract.receipt_path(self.work, f'prepared/{self.jobs[0][0]}/{self.jobs[0][2]}').unlink()
        with self.assertRaisesRegex(ValueError, 'Unreceipted'):
            contract.check_recovery_outputs(self.work)

    def test_finalization_failure_keeps_every_receipted_raw_estimate(self):
        self.run_reuse()
        for job in self.jobs:
            destination = self.work / job[4] / 'nnls_unweighted_fractions.tsv'
            destination.write_bytes((self.work / job[4] / 'epidish_rpc_fractions.tsv').read_bytes())
            contract.write_receipt(self.work, f'nnls/{job[0]}/{job[2]}', [f'{job[4]}/nnls_unweighted_fractions.tsv'])
        # Source input revalidation fails after all six comparison figures are written.
        result = subprocess.run(['Rscript', str(SCRIPTS / 'finalize_comparisons.R'),
                                 str(self.work), '2026-10-28_08-17-contract-test-v1'],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(list(self.work.glob('.finalization_attempts/1/comparisons/*/*/composition_comparison.pdf'))), 6, result.stderr)
        for job in self.jobs:
            self.assertTrue(contract.check_receipt(self.work, f'matched/{job[0]}/{job[2]}',
                                                   [f'{job[4]}/.uxm_weighted_matched.raw.csv']))
        self.assertFalse((self.work / 'EPIDISHR_COMPLETED').exists())
        self.assertEqual(list(self.work.glob('comparisons/*/*/composition_comparison.pdf')), [])
        contract.check_recovery_outputs(self.work)

    def test_published_finalization_outputs_without_receipt_block_recovery(self):
        put(self.work / 'comparisons/summary/marker_counts.tsv', 'partial\n')
        before = (self.work / 'provenance/execution_attempts.json').read_bytes()
        with patch.object(contract, 'build_contract', return_value=json.loads((self.work / 'provenance/recovery_contract.json').read_text())):
            with self.assertRaisesRegex(ValueError, 'Unreceipted partial stage finalization'):
                contract.reserve(Namespace(work_dir=self.work, recover=True,
                                           started_utc='2026-10-28T10:25:31Z', attempt_id='2'))
        self.assertEqual(before, (self.work / 'provenance/execution_attempts.json').read_bytes())





    def test_unreceipted_auxiliary_provenance_blocks_recovery_before_attempt_append(self):
        for name in ('unweighted_nnls.json', 'matched_uxm_reuse.json'):
            with self.subTest(name=name):
                path = self.work / 'provenance' / name
                write_json(path, {'partial': True})
                timing = (self.work / 'provenance/execution_attempts.json').read_bytes()
                with patch.object(contract, 'build_contract', return_value=json.loads((self.work / 'provenance/recovery_contract.json').read_text())):
                    with self.assertRaisesRegex(ValueError, 'Unreceipted partial stage'):
                        contract.reserve(Namespace(work_dir=self.work, recover=True,
                                                   started_utc='2026-10-28T10:25:31Z', attempt_id='2'))
                self.assertEqual(timing, (self.work / 'provenance/execution_attempts.json').read_bytes())
                self.assertEqual(json.loads(path.read_text()), {'partial': True})
                path.unlink()

    def test_standalone_plot_failure_records_start_and_end_without_completion(self):
        put(self.baseline / 'provenance/sample_metadata.tsv', 'sample\tpatient\tcompartment\tblast_pct\ns1\tp1\tBM cellular DNA\t50\ns2\tp1\tPB cellular DNA\t40\n')
        put(self.baseline / 'comparisons/summary/fractions_long.tsv.gz',
            'sample\tprofile\tcoverage_threshold\tmethod\tcell_type\tfraction\n'
            's1\tfull\t1\tepidish_rpc\tBlood-B\t0.7\n'
            's1\tfull\t1\tepidish_rpc\tBlood-T\t0.3\n'
            's2\tfull\t1\tepidish_rpc\tBlood-B\t0.2\n'
            's2\tfull\t1\tepidish_rpc\tBlood-T\t0.8\n')
        put(self.baseline / 'comparisons/summary/epidish_warnings_all.tsv', 'sample\twarning\tprofile\tcoverage_threshold\n')
        output = self.root / '2026-10-28_08-17-epidish-plot-failure-v1'
        # This failure-path fixture isolates the plotting error after reservation;
        # clean/dirty Git behavior has separate real temporary-repository tests.
        bindir = self.root / 'mock-git'; bindir.mkdir()
        mock_git = bindir / 'git'
        mock_git.write_text('#!/bin/sh\ncase "$*" in\n*status*) exit 0;;\n*--show-toplevel*) echo '+shlex.quote(str(SCRIPTS.parents[2]))+';;\n*) echo '+('a' * 40)+';;\nesac\n')
        mock_git.chmod(0o755)
        result = subprocess.run(['Rscript', str(SCRIPTS / 'plot_epidish_standalone.R'),
                                 str(self.baseline), str(output)], capture_output=True, text=True,
                                 env=dict(os.environ, PATH=str(bindir)+':'+os.environ['PATH']))
        self.assertNotEqual(result.returncode, 0)
        evidence = json.loads((output / 'epidish_plot_provenance.json').read_text())
        self.assertEqual(evidence['status'], 'failed')
        self.assertIn('started_utc', evidence)
        self.assertIn('ended_utc', evidence)
        self.assertNotIn('completed_utc', evidence)
        invalid = self.root / '2026-02-30_08-17-epidish-invalid-v1'
        result = subprocess.run(['Rscript', str(SCRIPTS / 'plot_epidish_standalone.R'),
                                 str(self.baseline), str(invalid)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(invalid.exists())

    def test_frozen_runtime_identity_is_exact_and_attempt_ids_are_distinct(self):
        checkout = self.root / 'runtime-uxm'
        checkout.mkdir()
        subprocess.run(['git', 'init', '-q', str(checkout)], check=True)
        subprocess.run(['git', '-C', str(checkout), '-c', 'user.name=Fixture',
                        '-c', 'user.email=fixture@example.invalid', 'commit', '--allow-empty', '-q', '-m', 'Fixture'], check=True)
        for label in ('full-runtime', 'heme-runtime'):
            run = self.root / label
            put(run / 'pat_manifest.tsv', 'sample\tpat\tcsi\tlbeta\ns1\tpat/s1.pat.gz\tpat/s1.pat.gz.csi\tlbeta/s1.lbeta\n')
            put(run / 'results/lineage_fractions.csv', 'CellType,s1\nBlood-B,1\n')
            for name in ('pat/s1.pat.gz', 'pat/s1.pat.gz.csi', 'lbeta/s1.lbeta'):
                put(run / name, 'fixture\n')
        put(self.root / 'full-runtime/intermediates/l4/s1.mem.test.homog.gz', 'fixture\n')
        for name in ('atlas-runtime.tsv', 'metadata-runtime.tsv', 'samplesheet-runtime.csv'):
            put(self.root / name, 'fixture\n')
        args = Namespace(work_dir=self.root / 'runtime-work', recover=False,
            atlas=self.root / 'atlas-runtime.tsv', metadata=self.root / 'metadata-runtime.tsv',
            samplesheet=self.root / 'samplesheet-runtime.csv', full_run=self.root / 'full-runtime',
            heme_run=self.root / 'heme-runtime', uxm_checkout=str(checkout), source_results=self.root,
            r_identity='4.5 EpiDISH 2.26.0', rpc_maxit=75, baseline='',
            started_utc='2026-10-28T08:17:21Z', attempt_id='101')
        with patch.object(contract, 'version', return_value='fixture-v1'):
            contract.reserve(args)
            frozen = (args.work_dir / 'provenance/recovery_contract.json').read_bytes()
            args.recover = True
            args.attempt_id = '102'
            args.started_utc = '2026-10-28T10:25:31Z'
            contract.reserve(args)
            self.assertEqual(frozen, (args.work_dir / 'provenance/recovery_contract.json').read_bytes())
            with self.assertRaisesRegex(ValueError, 'attempt ID already exists'):
                contract.reserve(args)
        before = (args.work_dir / 'provenance/execution_attempts.json').read_bytes()
        args.attempt_id = '103'
        with patch.object(contract, 'version', return_value='fixture-v2'):
            with self.assertRaisesRegex(ValueError, 'contract differs'):
                contract.reserve(args)
        self.assertEqual(before, (args.work_dir / 'provenance/execution_attempts.json').read_bytes())

    def test_complete_synthetic_finalization_publishes_receipted_sections(self):
        checkout = self.root / 'fixture-uxm'
        checkout.mkdir()
        subprocess.run(['git', 'init', '-q', str(checkout)], check=True)
        subprocess.run(['git', '-C', str(checkout), '-c', 'user.name=Fixture',
                        '-c', 'user.email=fixture@example.invalid', 'commit', '--allow-empty', '-q', '-m', 'Fixture'], check=True)
        commit = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
        inputs = copy.deepcopy(self.inputs)
        for label in ('atlas', 'environment'):
            path = self.root / f'{label}.txt'
            put(path, 'synthetic fixture\n')
            inputs[label] = contract.record(path)
        metadata = self.work / 'provenance/sample_metadata.tsv'
        put(metadata, 'sample\tpatient\tcompartment\tblast_pct\ns1\tp1\tBM cellular DNA\t50\ns2\tp1\tPB cellular DNA\t40\n')
        inputs['metadata'] = contract.record(metadata)
        inputs['uxm_checkout'] = {'path': str(checkout), 'commit': commit}
        inputs['r_version'] = subprocess.check_output(['Rscript', '-e', 'cat(as.character(getRversion()))'], text=True).strip()
        inputs['epidish_version'] = '2.26.0'
        for label in ('full_uxm_run', 'hematopoietic_uxm_run'):
            run = self.root / label
            put(run / 'pat_manifest.tsv', 'sample\tpat\tcsi\tlbeta\ns1\tp1\tc1\tl1\ns2\tp2\tc2\tl2\n')
            put(run / 'results/lineage_fractions.csv', (self.root / 'native.csv').read_text())
            inputs[label] = {'path': str(run), 'pat_manifest_sha256': contract.sha256_file(run / 'pat_manifest.tsv'),
                             'native_result_sha256': contract.sha256_file(run / 'results/lineage_fractions.csv')}
        inputs['regional_tables'] = []
        for sample in ('s1', 's2'):
            path = self.root / f'{sample}.homog.gz'
            put(path, 'synthetic regional table\n')
            inputs['regional_tables'].append({'sample': sample, **contract.record(path)})
        inputs['pat_artifacts'] = []
        for sample in ('s1', 's2'):
            for kind in ('pat', 'csi', 'lbeta'):
                path = self.root / f'{sample}.{kind}'
                put(path, 'synthetic artifact\n')
                inputs['pat_artifacts'].append({'sample': sample, 'kind': kind, **contract.record(path)})
        write_json(self.work / 'provenance/input_manifest.json', inputs)
        contract.receipt_path(self.work, 'preparation').unlink()
        contract.write_receipt(self.work, 'preparation', [
            'provenance/comparison_jobs.tsv', 'provenance/pat_manifest.tsv',
            'provenance/sample_metadata.tsv', 'provenance/input_manifest.json'])
        for job in self.jobs:
            raw = f'{job[4]}/.uxm_weighted_matched.raw.csv'
            put(self.work / raw, (self.root / 'native.csv').read_text())
            contract.write_receipt(self.work, f'matched/{job[0]}/{job[2]}', [raw])
            nnls = f'{job[4]}/nnls_unweighted_fractions.tsv'
            put(self.work / nnls, (self.work / job[4] / 'epidish_rpc_fractions.tsv').read_text())
            contract.write_receipt(self.work, f'nnls/{job[0]}/{job[2]}', [nnls])
        put(self.work / 'provenance/commands.log', 'synthetic finalization\n')
        write_json(self.work / 'provenance/unweighted_nnls.json', {'synthetic': True})
        contract.write_receipt(self.work, 'nnls_provenance', ['provenance/unweighted_nnls.json'])
        command = ['Rscript', str(SCRIPTS / 'finalize_comparisons.R'), str(self.work),
                   '2026-10-28_08-17-contract-test-v1']
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((self.work / 'EPIDISHR_COMPLETED').is_file())
        self.assertTrue(contract.check_receipt(self.work, 'finalization', contract.finalization_outputs()))
        marker = (self.work / 'EPIDISHR_COMPLETED').read_bytes()
        receipt = contract.receipt_path(self.work, 'finalization').read_bytes()
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(marker, (self.work / 'EPIDISHR_COMPLETED').read_bytes())
        self.assertEqual(receipt, contract.receipt_path(self.work, 'finalization').read_bytes())
        contract.check_recovery_outputs(self.work)
        # Simulate interruption after receipt publication but before the completion marker.
        (self.work / 'EPIDISHR_COMPLETED').unlink()
        original_commands = (self.work / 'provenance/commands.log').read_bytes()
        with patch.object(contract, 'build_contract', return_value=json.loads((self.work / 'provenance/recovery_contract.json').read_text())):
            contract.reserve(Namespace(work_dir=self.work, recover=True,
                                       started_utc='2026-10-28T10:25:31Z', attempt_id='2'))
        launcher = (SCRIPTS.parent / 'run_epidish_comparison.sh').read_text()
        start = launcher.index('# Completed finalization seals')
        selection = launcher[start:launcher.index('\nexec > >(tee', start)]
        subprocess.run(['bash', '-c', 'set -euo pipefail\n' +
            f'WORK_DIR={shlex.quote(str(self.work))}\nPROJECT_DIR={shlex.quote(str(self.root))}\nJOB_ID=2\n' + selection +
            '\nprintf "recovery attempt\\n" >> "${COMMANDS_LOG}"'], check=True)
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(original_commands, (self.work / 'provenance/commands.log').read_bytes())
        manifest = json.loads((self.work / 'provenance/run_manifest.json').read_text())
        self.assertEqual(manifest['commands_log_sha256'], contract.sha256_file(self.work / 'provenance/commands.log'))
        self.assertTrue((self.work / 'provenance/attempt_commands/2.log').is_file())
        contract.check_recovery_outputs(self.work)
        timing = (self.work / 'provenance/execution_attempts.json').read_bytes()
        put(self.work / 'provenance/commands.log', 'truncated\n')
        with patch.object(contract, 'build_contract', return_value=json.loads((self.work / 'provenance/recovery_contract.json').read_text())):
            with self.assertRaisesRegex(ValueError, 'provenance checksum mismatch: commands.log'):
                contract.reserve(Namespace(work_dir=self.work, recover=True,
                                           started_utc='2026-10-28T11:25:31Z', attempt_id='3'))
        self.assertEqual(timing, (self.work / 'provenance/execution_attempts.json').read_bytes())
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('provenance checksum mismatch', result.stderr)

    def publication_shell(self, source, final, between=''):
        launcher = (SCRIPTS.parent / 'run_epidish_comparison.sh').read_text()
        lock = launcher[launcher.index('command -v flock'):launcher.index('source "$(conda info --base)')]
        publication = launcher[launcher.index('"${CONTRACT_PYTHON}" "${OUTPUT_CONTRACT}" publish'):launcher.index('attempt_root="${FINAL_DIR}"')]
        assignments = {'WORK_DIR': source, 'WORK_BASE': source.parent, 'FINAL_DIR': final, 'RESULTS_BASE': final.parent,
                       'RECOVER_JOB_ID': '', 'CONTRACT_PYTHON': sys.executable,
                       'OUTPUT_CONTRACT': SCRIPTS.parents[1] / 'output_contract.py'}
        guard = SCRIPTS.parents[2] / 'scripts/lib/output_guard.py'
        reservation = " ".join(shlex.quote(str(x)) for x in (sys.executable, guard, "reserve", "--parent", final.parent, "--name", final.name, "--identity", json.dumps({"fixture": True}))) + " >/dev/null"
        return ('set -euo pipefail\n' + ''.join(f'{key}={shlex.quote(str(value))}\n' for key, value in assignments.items()) +
                lock + '\n' + reservation + '\n' + between + '\n' + publication)

    def test_concurrent_publications_share_destination_lock(self):
        final = self.root / 'published' / '2026-11-02_13-00-fixture-v1'
        first_work = self.root / 'work-first'
        second_work = self.root / 'work-second'
        first = subprocess.Popen(['bash', '-c', self.publication_shell(first_work, final,
            'mkdir "${WORK_DIR}"; echo READY; read -r release')],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.addCleanup(lambda: first.kill() if first.poll() is None else None)
        self.assertEqual(first.stdout.readline().strip(), 'READY')
        second = subprocess.run(['bash', '-c', self.publication_shell(second_work, final,
            'mkdir "${WORK_DIR}"')], capture_output=True, text=True)
        self.assertNotEqual(second.returncode, 0)
        self.assertIn('locked', second.stderr)
        self.assertFalse(second_work.exists())
        stdout, stderr = first.communicate('release\n', timeout=10)
        self.assertEqual(first.returncode, 0, stdout + stderr)
        self.assertTrue(final.is_dir())
        self.assertTrue(first_work.is_dir())
        self.assertEqual(list(first_work.iterdir()), [])

    def test_publication_destinations_and_late_collisions_preserve_estimates(self):
        self.run_reuse()
        final = self.root / 'published' / '2026-11-02_13-00-fixture-v1'
        for kind in ('file', 'directory', 'dangling'):
            for late in (False, True):
                with self.subTest(kind=kind, late=late):
                    target = final.with_name(f'2026-11-02_13-00-{kind}-{late}-v1')
                    target.parent.mkdir(parents=True, exist_ok=True)
                    collision = target / 'provenance' if late else target
                    if kind == 'file':
                        create = f'printf preserved > {shlex.quote(str(collision))}'
                    elif kind == 'directory':
                        create = f'mkdir {shlex.quote(str(collision))}'
                    else:
                        create = f'ln -s missing-target {shlex.quote(str(collision))}'
                    if not late:
                        subprocess.run(['bash', '-c', create], check=True)
                    result = subprocess.run(['bash', '-c', self.publication_shell(
                        self.work, target, create if late else '')], capture_output=True, text=True)
                    self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertTrue(os.path.lexists(collision))
                    if kind == 'file':
                        self.assertEqual(collision.read_text(), 'preserved')
                    elif kind == 'dangling':
                        self.assertEqual(os.readlink(collision), 'missing-target')
                    self.assertTrue(self.work.is_dir())
                    for job in self.jobs:
                        self.assertTrue(contract.check_receipt(self.work, f'matched/{job[0]}/{job[2]}',
                                                               [f'{job[4]}/.uxm_weighted_matched.raw.csv']))



if __name__ == '__main__':
    unittest.main()
