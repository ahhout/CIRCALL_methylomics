"""Safety tests use tiny disposable /tmp fixtures, never scientific files."""
import argparse
import contextlib
import importlib.util
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('cleanup', Path(__file__).resolve().parents[1] / 'scripts/cleanup_methylseq.py')
cleanup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cleanup)


class CleanupTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='emseq-cleanup-test-', dir='/tmp')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.repo = self.root / 'repo'; self.repo.mkdir()
        self.bulk = self.root / 'results/bulk'; self.bulk.mkdir(parents=True)
        self.run = self.bulk / sorted(cleanup.JULY)[0]; self.run.mkdir()
        for key, value in {'REPO': self.repo, 'BULK': self.bulk, 'RESULTS': self.bulk.parent,
                           'SCAN_ROOTS': (self.repo,), 'SCRATCH_METADATA': Path('/tmp')}.items():
            manager = patch.object(cleanup, key, value); manager.start(); self.addCleanup(manager.stop)
        manager = patch.object(cleanup, 'scheduler', return_value={'checked_utc': cleanup.now(), 'jobs': [], 'error': None})
        self.scheduler = manager.start(); self.addCleanup(manager.stop)
        self.make('bismark/deduplicated/S1.deduplicated.sorted.bam', b'bam-data')
        self.make('bismark/deduplicated/S1.deduplicated.sorted.bam.bai', b'index')
        self.make('bismark/methylation_calls/methylation_calls/CpG_OT_S1.txt.gz', b'calls')
        for name in ('pipeline_info/params.json', 'bismark/methylation_calls/methylation_calls/S1_splitting_report.txt',
                     'bismark/methylation_calls/methylation_calls/CpG_S1_splitting_report.txt.gz',
                     'bismark/methylation_calls/methylation_coverage/S1.cov.gz',
                     'bismark/coverage2cytosine/reports/S1.CpG_report.txt.gz',
                     'bismark/methylation_calls/methylation_calls/S1.M-bias.txt',
                     'raw/S1.fastq.gz', 'pat/S1.pat.gz', 'references/genome.fa',
                     '.nextflow/task/.command.sh', 'unknown.bin'):
            self.make(name, b'preserve')
        self.external = self.root / 'external-provenance.txt'; self.external.write_text('kept log')
        self.review_file = self.root / 'review.json'
        self.write_review()
        self.output = self.root / 'plan'
        self.args = argparse.Namespace(run_dir=[str(self.run)], artifacts='bam', manifest_dir=str(self.output),
                    dependency_review=str(self.review_file), apply=False, approved_manifest=None,
                    confirm_manifest_sha256=None)

    def make(self, name, content):
        path = self.run / name; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(content)
        return path

    def write_review(self, findings=None, jobs=None, run=None):
        run = run or self.run
        review = {'checked_utc': cleanup.now(), 'runs': {str(run): {
            'no_active_consumers': True, 'evidence': 'Fixture has no scientific consumers',
            'recovery_policy': 'accept_no_guaranteed_replay',
            'historical_metadata': {key: 'fixture' for key in ('pipeline_version', 'params', 'config', 'samplesheet',
                                                            'reference', 'scientific_settings', 'resume', 'missing_originals')},
            'external_records': [{'path': str(self.external), 'sha256': cleanup.hash_file(self.external)}]}},
            'findings': findings or {}, 'unrelated_jobs': jobs or {}}
        self.review_file.write_text(json.dumps(review))
        return review

    def plan(self):
        with contextlib.redirect_stdout(io.StringIO()): cleanup.make_plan(self.args)
        self.args.approved_manifest = str(self.output / 'plan.tsv')
        self.args.confirm_manifest_sha256 = cleanup.hash_file(self.output / 'plan.tsv')
        return json.loads((self.output / 'summary.json').read_text())

    def apply(self):
        with contextlib.redirect_stdout(io.StringIO()): cleanup.apply_plan(self.args)

    def snapshot(self):
        return {str(p.relative_to(self.run)): (p.read_bytes() if not p.is_symlink() else os.readlink(p))
                for p in cleanup.walk_files(self.run)}

    def test_default_dry_run_preserves_every_source_byte(self):
        before = self.snapshot(); summary = self.plan()
        self.assertEqual(self.snapshot(), before)
        self.assertEqual(summary['selected']['files'], 2)
        self.assertEqual(summary['mode'], 'dry-run')
        self.assertFalse((self.output / 'deletions.jsonl').exists())

    def test_refuse_broad_relative_unknown_and_protected_runs(self):
        for name in ('2026-08-31-emseq-integrated-8samples-trim-8-3-3-10-v1', 'work_emseq_integrated', 'unknown'):
            (self.bulk / name).mkdir()
        cases = [str(self.bulk), str(self.bulk.parent), 'relative', str(self.run / 'bismark'),
                 str(self.run) + '/../' + self.run.name, str(self.bulk / '*')]
        cases += [str(p) for p in self.bulk.iterdir() if p != self.run]
        for value in cases:
            with self.subTest(path=value), self.assertRaises(cleanup.CleanupError): cleanup.validate_runs([value])

    def test_e20_e21_require_migration_closure_and_file_specific_release(self):
        gates = ('e31_downstream_migration_verified', 'dmrichr_local_closure_verified',
                 'retained_consumers_independent_of_candidates')
        for name in sorted(cleanup.E20_E21):
            run = self.bulk / name
            run.mkdir()
            self.assertEqual(cleanup.validate_runs([str(run)]), [run])
            review = self.write_review(run=run)
            queue = {'error': None, 'jobs': []}
            self.assertEqual(len(cleanup.review_blockers(review, [run], [], queue)), 3)
            for gate in gates:
                review['runs'][str(run)][gate] = True
            self.assertEqual(cleanup.review_blockers(review, [run], [], queue), [])
            for gate in gates:
                review['runs'][str(run)][gate] = False
                reasons = cleanup.review_blockers(review, [run], [], queue)
                self.assertEqual(len(reasons), 1)
                self.assertIn(gate, reasons[0])
                review['runs'][str(run)][gate] = True

    def test_alias_and_symlinked_parent_refused(self):
        alias = self.root / 'alias'; alias.symlink_to(self.run, target_is_directory=True)
        parent_alias = self.root / 'bulk-alias'; parent_alias.symlink_to(self.bulk, target_is_directory=True)
        for path in (alias, parent_alias / self.run.name):
            with self.assertRaises(cleanup.CleanupError): cleanup.validate_runs([str(path)])

    def test_directory_symlink_never_traversed(self):
        outside = self.root / 'outside'; outside.mkdir(); (outside / 'keep.bam').write_text('external')
        (self.run / 'external').symlink_to(outside, target_is_directory=True)
        before = self.snapshot(); self.plan(); self.apply()
        self.assertTrue((self.run / 'external').is_symlink())
        self.assertEqual((outside / 'keep.bam').read_text(), 'external')

    def test_bam_symlink_and_orphan_bam_preserved(self):
        path = self.make('bismark/alignments/orphan.bam', b'orphan')
        alias = self.run / 'bismark/alignments/linked.bam'; alias.symlink_to(self.external)
        self.make('bismark/alignments/linked.bam.bai', b'index')
        summary = self.plan(); self.apply()
        self.assertEqual(summary['selected']['files'], 2)
        self.assertTrue(path.exists()); self.assertTrue(alias.is_symlink())
        self.assertTrue(Path(str(alias)+'.bai').exists())

    def test_hardlinked_candidate_refused(self):
        os.link(self.run / 'bismark/deduplicated/S1.deduplicated.sorted.bam', self.root / 'hardlink')
        with self.assertRaisesRegex(cleanup.CleanupError, 'multiple hard links'): self.plan()

    def test_duplicate_run_arguments_not_double_counted(self):
        self.args.run_dir *= 2
        summary = self.plan(); self.assertEqual(summary['selected']['files'], 2)
        self.assertEqual(summary['selected']['unique_inodes'], 2)

    def test_only_selected_classes_are_removed_and_journaled(self):
        self.args.artifacts = 'bam,context-calls'
        before = self.snapshot(); summary = self.plan(); self.apply()
        self.assertEqual(summary['selected']['files'], 3)
        rows = [json.loads(line) for line in (self.output / 'deletions.jsonl').read_text().splitlines()]
        self.assertEqual([r['event'] for r in rows], ['intent', 'unlinked'] * 3)
        removed = {r['path'] for r in rows if r['event'] == 'unlinked'}
        for name, contents in before.items():
            if name not in removed: self.assertEqual((self.run / name).read_bytes(), contents)
        self.assertTrue((self.run / 'PAYLOAD_RETENTION.md').is_file())
        result = json.loads((self.output / 'execution_summary.json').read_text())
        self.assertTrue(result['retained_verified'])
        self.assertEqual(result['removed']['allocated_bytes'], summary['selected']['allocated_bytes'])

    def test_report_and_reference_like_bam_names_retained(self):
        for name in ('reference.bam', 'genome.bam', 'cache.bam'):
            self.make('bismark/alignments/'+name, b'keep')
            self.make('bismark/alignments/'+name+'.bai', b'keep')
        summary = self.plan(); self.assertEqual(summary['selected']['files'], 2)

    def test_modified_plan_hash_refused(self):
        self.plan(); before = self.snapshot()
        with (self.output / 'plan.tsv').open('a') as handle: handle.write('changed\n')
        with self.assertRaisesRegex(cleanup.CleanupError, 'manifest hash mismatch'): self.apply()
        self.assertEqual(before, self.snapshot())

    def test_changed_retention_context_is_hash_bound(self):
        self.plan()
        context = self.output / 'context.json'
        context.write_text(context.read_text() + ' ')
        with self.assertRaisesRegex(cleanup.CleanupError, 'context changed'): self.apply()

    def test_changed_retained_checksum_companion_refused(self):
        self.plan(); (self.output / 'retained.sha256').write_text('wrong')
        with self.assertRaisesRegex(cleanup.CleanupError, 'companion record changed'): self.apply()

    def test_changed_selected_file_refused(self):
        self.plan(); self.make('bismark/deduplicated/S1.deduplicated.sorted.bam', b'new much larger BAM')
        with self.assertRaisesRegex(cleanup.CleanupError, 'contents changed'): self.apply()
        self.assertFalse((self.output / 'deletions.jsonl').exists())

    def test_changed_retained_file_same_size_and_timestamp_refused(self):
        self.plan(); p = self.run / 'pipeline_info/params.json'; old = p.stat()
        p.write_bytes(b'changed!'); os.utime(p, ns=(old.st_atime_ns, old.st_mtime_ns))
        with self.assertRaisesRegex(cleanup.CleanupError, 'contents changed'): self.apply()

    def test_unreviewed_dependency_blocks_apply(self):
        (self.repo / 'consumer.conf').write_text(str(self.run))
        summary = self.plan(); self.assertTrue(summary['blocked'])
        with self.assertRaisesRegex(cleanup.CleanupError, 'unresolved blockers'): self.apply()

    def test_new_consumer_after_approval_blocks_apply(self):
        self.plan(); (self.repo / 'consumer.conf').write_text(str(self.run))
        with self.assertRaisesRegex(cleanup.CleanupError, 'Dependency evidence changed'): self.apply()

    def test_new_scheduler_job_blocks_apply(self):
        self.plan(); self.scheduler.return_value = {'checked_utc': cleanup.now(), 'error': None,
            'jobs': [{'id': '123', 'sha256': 'abc', 'description': 'new job'}]}
        with self.assertRaisesRegex(cleanup.CleanupError, 'Unreviewed current Slurm job'): self.apply()

    def test_failed_scheduler_check_blocks_apply(self):
        self.plan(); self.scheduler.return_value = {'checked_utc': cleanup.now(), 'jobs': [], 'error': 'unavailable'}
        with self.assertRaisesRegex(cleanup.CleanupError, 'scheduler state'): self.apply()

    def test_review_expiry_requires_new_manifest(self):
        review = self.write_review(); review['checked_utc'] = '2000-01-01T00:00:00+00:00'
        self.review_file.write_text(json.dumps(review)); self.plan()
        with self.assertRaisesRegex(cleanup.CleanupError, 'unresolved blockers'): self.apply()

    def test_missing_review_is_dry_run_only(self):
        self.args.dependency_review = None
        self.assertTrue(self.plan()['blocked'])
        with self.assertRaisesRegex(cleanup.CleanupError, 'unresolved blockers'): self.apply()

    def test_e27_requires_explicit_mbias_closure(self):
        other = self.bulk / cleanup.E27; self.run.rename(other); self.run = other
        self.args.run_dir = [str(other)]; self.write_review(run=other)
        self.assertTrue(any('M-bias' in reason for reason in self.plan()['blocked']))
        with self.assertRaises(cleanup.CleanupError): self.apply()

    def test_no_plan_overwrite_or_implicit_retry(self):
        self.plan()
        with self.assertRaisesRegex(cleanup.CleanupError, 'overwrite'): self.plan()
        self.apply()
        with self.assertRaisesRegex(cleanup.CleanupError, 'attempt already exists'): self.apply()

    def test_scope_or_class_change_refused(self):
        self.plan(); self.args.artifacts = 'bam,context-calls'
        with self.assertRaisesRegex(cleanup.CleanupError, 'scope'): self.apply()

    def test_interruption_is_journaled_and_cannot_auto_resume(self):
        self.plan()
        real_unlink = cleanup.unlink_checked
        calls = []
        def interrupted(row):
            calls.append(row)
            if len(calls) == 2: raise KeyboardInterrupt('fixture interruption')
            real_unlink(row)
        with patch.object(cleanup, 'unlink_checked', side_effect=interrupted):
            with self.assertRaisesRegex(cleanup.CleanupError, 'Apply stopped'): self.apply()
        events = [json.loads(line)['event'] for line in (self.output / 'deletions.jsonl').read_text().splitlines()]
        self.assertEqual(events, ['intent', 'unlinked', 'intent', 'failed'])
        result = json.loads((self.output / 'execution_summary.json').read_text())
        self.assertEqual(result['removed']['files'], 1); self.assertEqual(result['status'], 'failed')
        with self.assertRaisesRegex(cleanup.CleanupError, 'attempt already exists'): self.apply()

    def test_selected_parent_symlink_swap_refused_at_unlink(self):
        candidates, _, _ = cleanup.inventory([self.run], ['bam'])
        parent = self.run / 'bismark/deduplicated'; moved = self.root / 'moved'; parent.rename(moved)
        parent.symlink_to(moved, target_is_directory=True)
        with self.assertRaises(OSError): cleanup.unlink_checked(candidates[0])
        self.assertTrue((moved / Path(candidates[0]['path']).name).exists())

    def test_changed_external_provenance_blocks_apply(self):
        self.plan(); self.external.write_text('changed external log')
        with self.assertRaisesRegex(cleanup.CleanupError, 'External retained evidence changed'): self.apply()

    def test_missing_pipeline_provenance_blocks_apply(self):
        (self.run / 'pipeline_info/params.json').unlink()
        self.assertTrue(any('pipeline_info' in reason for reason in self.plan()['blocked']))
        with self.assertRaisesRegex(cleanup.CleanupError, 'unresolved blockers'): self.apply()

    def test_explicit_review_clears_historical_reference_and_unrelated_job(self):
        (self.repo / 'historical.md').write_text('Historical input: '+str(self.run))
        findings = cleanup.scan_dependencies([self.run])
        job = {'id': '123', 'sha256': 'fingerprint', 'description': 'unrelated fixture allocation'}
        self.scheduler.return_value = {'checked_utc': cleanup.now(), 'error': None, 'jobs': [job]}
        self.write_review(findings={f['id']: {'status': 'historical', 'reason': 'Retained provenance only'} for f in findings},
                          jobs={'123': {'sha256': job['sha256'], 'reason': 'Independent fixture job'}})
        self.assertEqual(self.plan()['blocked'], [])
        self.apply()

    def test_directory_inspection_error_does_not_silently_skip(self):
        with self.assertRaisesRegex(cleanup.CleanupError, 'Cannot inspect directory'):
            cleanup.walk_error(PermissionError('fixture unreadable directory'))

    def test_sync_failure_after_unlink_does_not_underreport_removed_file(self):
        self.plan()
        with patch.object(cleanup, 'sync_unlink_parent', side_effect=OSError('fixture sync failure')):
            with self.assertRaisesRegex(cleanup.CleanupError, 'Apply stopped'): self.apply()
        events = [json.loads(line)['event'] for line in (self.output / 'deletions.jsonl').read_text().splitlines()]
        self.assertEqual(events, ['intent', 'unlinked', 'durability-error'])
        result = json.loads((self.output / 'execution_summary.json').read_text())
        self.assertEqual(result['removed']['files'], 1)
        self.assertEqual(result['status'], 'failed')

    def test_manifest_output_cannot_pollute_results_or_work_cache(self):
        for path in (self.run / 'manifest', self.bulk.parent / 'manifest',
                     cleanup.SCRATCH_METADATA / 'em_seq/work_emseq_integrated/manifest'):
            with self.assertRaises(cleanup.CleanupError): cleanup.validate_metadata_destination(path)


if __name__ == '__main__':
    unittest.main()
