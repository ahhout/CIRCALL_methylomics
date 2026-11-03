import concurrent.futures
import errno
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('output_guard', Path(__file__).resolve().parents[1] / 'scripts/lib/output_guard.py')
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)
NAME = '2026-11-02_14-30-test-v1'
IDENTITY = {'project_commit': 'a' * 40, 'pipeline_commit': 'b' * 40, 'input_sha256': 'c' * 64}


class OutputGuardTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir='/tmp')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_reject_every_existing_kind_without_changes(self):
        for kind in ['empty', 'nonempty', 'file', 'symlink', 'dangling']:
            parent = self.root / kind; parent.mkdir()
            dest = parent / NAME
            if kind in ['empty', 'nonempty']:
                dest.mkdir()
                if kind == 'nonempty': (dest / 'sentinel').write_text('original')
            elif kind == 'file': dest.write_text('original')
            elif kind == 'symlink': dest.symlink_to(self.root, target_is_directory=True)
            else: dest.symlink_to(parent / 'missing')
            before = dest.lstat()
            with self.assertRaises(FileExistsError): guard.reserve(parent, NAME, IDENTITY)
            self.assertEqual(before.st_ino, dest.lstat().st_ino)
            if kind == 'nonempty': self.assertEqual((dest / 'sentinel').read_text(), 'original')

    def test_concurrent_reservations_have_one_winner(self):
        def attempt(_):
            try: guard.reserve(self.root, NAME, IDENTITY); return True
            except FileExistsError: return False
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            self.assertEqual(sum(pool.map(attempt, range(16))), 1)
        with self.assertRaises(FileExistsError): guard.reserve(self.root, NAME, IDENTITY)

    def test_resume_exact_identity_only_and_completed_immutable(self):
        dest = guard.reserve(self.root, NAME, IDENTITY)
        with self.assertRaises(ValueError): guard.resume(dest, IDENTITY)
        guard.update(dest, status='running', session_id='session-a')
        guard.finish(dest, False)
        self.assertEqual(guard.resume(dest, IDENTITY, 'session-a')['status'], 'failed')
        with self.assertRaises(ValueError): guard.resume(dest, {**IDENTITY, 'input_sha256': 'changed'})
        with self.assertRaises(ValueError): guard.resume(dest, IDENTITY, 'session-b')
        with self.assertRaises(ValueError): guard.update(dest, session_id='session-b')
        guard.finish(dest, True)
        original = (dest / guard.RECORD).read_bytes()
        with self.assertRaises(ValueError): guard.resume(dest, IDENTITY)
        with self.assertRaises(ValueError): guard.update(dest, status='running')
        with self.assertRaises(FileExistsError): guard.reserve(self.root, NAME, IDENTITY)
        self.assertEqual((dest / guard.RECORD).read_bytes(), original)

    def test_failed_initialization_retains_reservation(self):
        with patch.object(guard, '_write', side_effect=OSError('fixture failure')):
            with self.assertRaises(OSError): guard.reserve(self.root, NAME, IDENTITY)
        self.assertTrue((self.root / NAME).is_dir())
        with self.assertRaises(FileNotFoundError): guard.resume(self.root / NAME, IDENTITY)
        with self.assertRaises(FileExistsError): guard.reserve(self.root, NAME, IDENTITY)

    def test_parent_alias_traversal_and_invalid_dates(self):
        alias = self.root / 'alias'; alias.symlink_to(self.root)
        with self.assertRaises(ValueError): guard.reserve(alias, NAME, IDENTITY)
        for name in ['../'+NAME, '2026-02-30_12-00-test-v1', '2026-11-02_25-00-test-v1', 'plain']:
            with self.assertRaises(ValueError): guard.reserve(self.root, name, IDENTITY)

    def test_links_never_repoint(self):
        a = guard.reserve(self.root, NAME, IDENTITY)
        b = guard.reserve(self.root, NAME+'-other', IDENTITY)
        link = self.root / 'navigation'
        guard.publish_link(link, a); guard.publish_link(link, a)
        with self.assertRaises(ValueError): guard.publish_link(link, b)
        self.assertEqual(link.resolve(), a)
        normal = self.root / 'normal'; normal.mkdir()
        with self.assertRaises(ValueError): guard.publish_link(normal, a)

    def test_cannot_create_nested_bundle_in_completed_output(self):
        output = guard.reserve(self.root, NAME, IDENTITY)
        nested = output / 'figures'; nested.mkdir()
        guard.finish(output, True)
        with self.assertRaises(ValueError): guard.reserve(nested, NAME+'-figure', IDENTITY)
        # Read/inspect the completed bundle itself remains allowed.
        self.assertEqual(guard.read_record(output)['status'], 'completed')

    def test_publish_contents_no_replace_and_completed_refused(self):
        output = guard.reserve(self.root, NAME, IDENTITY)
        stage = self.root / 'stage'; stage.mkdir(); (stage / 'result').write_text('new')
        (output / 'result').write_text('old')
        with self.assertRaises(FileExistsError): guard.publish_contents(stage, output)
        self.assertEqual((output / 'result').read_text(), 'old')
        (stage / 'result').rename(stage / 'other')
        guard.publish_contents(stage, output)
        self.assertEqual((output / 'other').read_text(), 'new')
        self.assertTrue(stage.is_dir())
        guard.finish(output, True)
        with self.assertRaises(ValueError): guard.publish_contents(stage, output)

    def test_beegfs_fallback_exclusive_hardlink_transfer(self):
        output = guard.reserve(self.root, NAME, IDENTITY)
        stage = self.root / 'stage'; stage.mkdir()
        nested = stage / 'data'; nested.mkdir(); (nested / 'file').write_text('payload')
        inode = (nested / 'file').stat().st_ino
        (stage / 'link').symlink_to('data/file')
        with patch.object(guard, '_rename_noreplace', side_effect=OSError(errno.EINVAL, 'fixture')):
            guard.publish_contents(stage, output)
        self.assertEqual((output / 'data/file').stat().st_ino, inode)
        self.assertEqual((output / 'link').read_text(), 'payload')
        self.assertTrue(nested.is_dir())
        self.assertFalse((nested / 'file').exists())
        (stage / 'link').write_text('conflict')
        with self.assertRaises(FileExistsError): guard.publish_contents(stage, output)
        self.assertEqual((output / 'data/file').read_text(), 'payload')

    def test_publication_cannot_remove_completed_source_or_its_children(self):
        previous = guard.reserve(self.root, NAME, IDENTITY)
        summary = previous / 'summary'; summary.mkdir()
        payload = summary / 'result'; payload.write_text('published')
        guard.finish(previous, True)
        output = guard.reserve(self.root, NAME+'-new', IDENTITY)
        original = (previous / guard.RECORD).read_bytes()
        for source in (previous, summary):
            with self.assertRaises(ValueError): guard.publish_contents(source, output)
        self.assertEqual(payload.read_text(), 'published')
        self.assertEqual((previous / guard.RECORD).read_bytes(), original)
        self.assertEqual([p.name for p in output.iterdir()], [guard.RECORD])

    def test_publication_rejects_source_alias_and_wrong_staging_receipt(self):
        output = guard.reserve(self.root, NAME, IDENTITY)
        stage = self.root / ('.'+NAME+'.running-1'); stage.mkdir()
        (stage / 'result').write_text('new')
        alias = self.root / 'alias'; alias.symlink_to(self.root)
        with self.assertRaises(ValueError): guard.publish_contents(alias / stage.name, output)
        receipt = stage / 'staging-ownership.json'
        receipt.write_text(json.dumps({'section': True, 'owner_output': 'different'}))
        with self.assertRaises(ValueError): guard.publish_contents(stage, output)
        self.assertTrue((stage / 'result').is_file())
        receipt.write_text(json.dumps({'section': True, 'owner_output': str(output)}))
        guard.publish_contents(stage, output)
        self.assertEqual((output / 'result').read_text(), 'new')

    def test_historical_completed_source_is_untouched(self):
        output = guard.reserve(self.root, NAME, IDENTITY)
        historical = self.root / 'historical'; historical.mkdir()
        (historical / 'COMPLETED').write_text('old completion')
        (historical / 'result').write_text('old result')
        with self.assertRaises(ValueError): guard.publish_contents(historical, output)
        self.assertEqual((historical / 'result').read_text(), 'old result')
        self.assertEqual((historical / 'COMPLETED').read_text(), 'old completion')
        self.assertEqual([p.name for p in output.iterdir()], [guard.RECORD])


if __name__ == '__main__': unittest.main()
