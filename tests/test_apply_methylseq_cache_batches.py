import copy
import json
from pathlib import Path
import sys
import tarfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import apply_methylseq_cache_batches as group
from test_cleanup_methylseq_cache import CacheTests, cache


class GroupTests(unittest.TestCase):
    def setUp(self):
        self.fixture = CacheTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.patcher = patch.object(group, 'c', cache)
        self.patcher.start(); self.addCleanup(self.patcher.stop)
        self.protected = self.fixture.make_task(self.fixture.roots[0], 'bb', '2' * 30)
        prefix, digest = sorted(cache.ABANDONED)[0].split('/')
        self.abandoned = self.fixture.make_task(self.fixture.roots[0], prefix, digest)
        (self.abandoned / '.exitcode').write_text('143')
        first, self.state = self.fixture.manifest()
        m = json.loads(first.read_text())
        m['batch'] = 'abandoned_alignment'
        m['tasks'] = [str(self.abandoned)]
        m['entries'] = self.state['tasks'][str(self.abandoned)]
        m['review']['tasks'] = {str(self.abandoned): {'no_active_consumers': True, 'evidence': 'fixture'}}
        second = first.parent / 'second.json'; second.write_bytes(cache.encode(m))
        self.paths = [str(first), str(second)]
        self.hashes = [cache.sha(Path(p).read_bytes()) for p in self.paths]
        self.output = self.fixture.base / 'group'

    def run_group(self, apply=True):
        with patch.object(group, 'collect_parallel', side_effect=self.fixture.state), \
             patch.object(cache, 'scheduler', return_value={'jobs': [], 'error': None}), \
             patch.object(cache, 'metadata_destination', side_effect=lambda p: Path(p)):
            group.run(self.paths, self.hashes, str(self.output), apply)

    def test_group_apply_and_verified_archive(self):
        self.run_group()
        result = json.loads((self.output / 'execution.json').read_text())
        self.assertEqual(result['tasks'], 2)
        self.assertFalse(self.fixture.task.exists()); self.assertFalse(self.abandoned.exists())
        self.assertEqual((self.protected / 'payload.bam').read_bytes(), b'payload')
        with tarfile.open(self.output / 'task_provenance.tar') as bundle:
            self.assertEqual(len(bundle.getmembers()), 3)
        events = [json.loads(x) for x in (self.output / 'deletions.jsonl').read_text().splitlines()]
        self.assertEqual(events[-1]['event'], 'complete')
        self.assertEqual(sum(e['event'] == 'intent' for e in events), result['removed_entries'])

    def test_dry_run_does_not_delete(self):
        self.run_group(False)
        self.assertTrue(self.fixture.task.exists())
        self.assertFalse((self.output / 'deletions.jsonl').exists())

    def test_changed_manifest_refused(self):
        Path(self.paths[0]).write_text(Path(self.paths[0]).read_text() + '\n')
        with self.assertRaisesRegex(cache.Refused, 'hash mismatch'):
            self.run_group()
        self.assertTrue(self.fixture.task.exists())

    def test_duplicate_batch_refused(self):
        with self.assertRaises(cache.Refused):
            group.load_manifests(self.paths * 2, self.hashes * 2)

    def test_added_file_refused(self):
        (self.fixture.task / 'extra').write_text('keep')
        with self.assertRaisesRegex(cache.Refused, 'contents changed'):
            self.run_group()
        self.assertTrue((self.fixture.task / 'payload.bam').exists())

    def test_interrupted_task_keeps_durable_intents(self):
        original = cache.remove_entry
        calls = []
        def fail(row):
            calls.append(row)
            if len(calls) == 2:
                raise OSError('injected interruption')
            original(row)
        with patch.object(cache, 'remove_entry', side_effect=fail):
            with self.assertRaisesRegex(OSError, 'injected'):
                self.run_group()
        events = [json.loads(x) for x in (self.output / 'deletions.jsonl').read_text().splitlines()]
        self.assertEqual(events[-1]['event'], 'stopped')
        self.assertEqual(sum(e['event'] == 'removed' for e in events), 1)
        self.assertGreater(sum(e['event'] == 'intent' for e in events), 1)
        self.assertFalse((self.output / 'execution.json').exists())


if __name__ == '__main__':
    unittest.main()
