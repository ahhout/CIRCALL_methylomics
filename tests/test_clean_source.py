"""Real temporary Git fixtures for clean committed-source execution gates."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts/lib'))
import output_guard


class CleanSource(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir='/tmp')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.git('init', '-q')
        self.git('config', 'user.name', 'Fixture')
        self.git('config', 'user.email', 'fixture@example.invalid')
        (self.root / 'source.py').write_text('x = 1\n')
        (self.root / '.gitignore').write_text('results/\n')
        self.git('add', 'source.py', '.gitignore')
        self.git('commit', '-qm', 'Fixture')

    def git(self, *args):
        return subprocess.check_output(['git', '-C', str(self.root), *args], text=True).strip()

    def test_clean_committed_file_and_ignored_outputs(self):
        revision = self.git('rev-parse', 'HEAD')
        (self.root / 'results').mkdir()
        (self.root / 'results/fixture').write_text('ignored generated output')
        self.assertEqual(output_guard.require_clean_source(self.root / 'source.py'), revision)

    def test_modified_source_rejected(self):
        (self.root / 'source.py').write_text('x = 2\n')
        with self.assertRaisesRegex(ValueError, 'clean committed source'):
            output_guard.require_clean_source(self.root)

    def test_staged_source_rejected(self):
        (self.root / 'source.py').write_text('x = 2\n')
        self.git('add', 'source.py')
        with self.assertRaisesRegex(ValueError, 'clean committed source'):
            output_guard.require_clean_source(self.root)

    def test_nonignored_untracked_source_rejected(self):
        (self.root / 'new-helper.py').write_text('x = 2\n')
        with self.assertRaisesRegex(ValueError, 'clean committed source'):
            output_guard.require_clean_source(self.root)

    def test_nonrepository_source_rejected(self):
        with tempfile.TemporaryDirectory(dir='/tmp') as directory:
            with self.assertRaisesRegex(ValueError, 'committed source in a Git repository'):
                output_guard.require_clean_source(directory)


if __name__ == '__main__':
    unittest.main()
