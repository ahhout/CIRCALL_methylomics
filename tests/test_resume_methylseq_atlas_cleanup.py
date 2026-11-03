import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import resume_methylseq_atlas_cleanup as recovery


class RecoveryTests(unittest.TestCase):
    def test_mount_mapping_only_changes_expected_device(self):
        row = {'path': '/path/to/hpc-storage/a', 'device': 42, 'inode': 99}
        normalized = recovery.mapped_devices(row, [('/path/to/hpc-storage', 43, 42)])
        self.assertEqual(normalized, dict(row, device=43))
        self.assertEqual(recovery.mapped_devices(normalized, [('/path/to/hpc-storage', 43, 42)], False), row)
        with self.assertRaises(recovery.cache.Refused):
            recovery.mapped_devices(dict(row, device=999), [('/path/to/hpc-storage', 43, 42)])

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.task = self.root / 'task'
        self.task.mkdir()
        self.payload = self.task / 'payload'
        self.payload.write_text('retired data')
        self.manifest = {'tasks': [str(self.task)], 'entries': [
            {'path': str(p), **recovery.cache.identity(p)} for p in (self.task, self.payload)]}
        self.journal = self.root / 'journal'
        self.events = [{'event': 'start'}]

    def check(self):
        self.journal.write_text(''.join(json.dumps(r) + '\n' for r in self.events))
        return recovery.reconcile(self.manifest, self.journal)

    def test_intent_without_unlink_remains(self):
        self.events.append({'event': 'intent', 'path': str(self.payload)})
        remaining, recovered = self.check()
        self.assertEqual(len(remaining), 2)
        self.assertEqual(recovered, [])

    def test_unlink_without_completion_is_reconciled(self):
        self.events.append({'event': 'intent', 'path': str(self.payload)})
        self.payload.unlink()
        remaining, recovered = self.check()
        self.assertEqual(len(remaining), 1)
        self.assertEqual(recovered, [str(self.payload)])

    def test_unrecorded_absence_refused(self):
        self.payload.unlink()
        with self.assertRaisesRegex(recovery.cache.Refused, 'without prior deletion intent'):
            self.check()

    def test_reappeared_path_refused(self):
        self.events += [{'event': k, 'path': str(self.payload)} for k in ('intent', 'removed')]
        with self.assertRaisesRegex(recovery.cache.Refused, 'reappeared'):
            self.check()

    def test_added_file_refused(self):
        (self.task / 'unexpected').touch()
        with self.assertRaisesRegex(recovery.cache.Refused, 'Unexpected directory contents'):
            self.check()

    def test_changed_payload_refused(self):
        self.payload.write_text('new content')
        with self.assertRaisesRegex(recovery.cache.Refused, 'changed'):
            self.check()

    def test_out_of_scope_journal_refused(self):
        self.events.append({'event': 'intent', 'path': str(self.root / 'other')})
        with self.assertRaisesRegex(recovery.cache.Refused, 'outside approved scope'):
            self.check()

    def test_symlink_directory_not_followed(self):
        linked = self.task / 'link'
        linked.symlink_to(self.root, target_is_directory=True)
        self.manifest['entries'].append({'path': str(linked), **recovery.cache.identity(linked)})
        remaining, recovered = self.check()
        self.assertEqual(len(remaining), 3)


if __name__ == '__main__':
    unittest.main()
