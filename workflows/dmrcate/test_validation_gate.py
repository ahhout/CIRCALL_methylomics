"""Exercise the launcher source gate with tiny inputs before any submission."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

WORKFLOW = Path(__file__).resolve().parent


class RuntimeValidationGate(unittest.TestCase):
    def test_stale_source_validation_fails_before_reservation_or_submission(self):
        launcher = (WORKFLOW / 'run_dmrcate.sh').read_text()
        prefix = launcher[:launcher.index('root="$(python -c')]
        with tempfile.TemporaryDirectory(dir='/tmp') as directory:
            root = Path(directory)
            validation = root / 'validation'
            (validation / 'r').mkdir(parents=True)
            runtime = root / 'runtime'; runtime.mkdir()
            (validation / 'VALIDATED').write_text('fixture')
            (validation / 'sif.sha256').write_text('a' * 64 + ' fixture.sif\n')
            (runtime / 'sif.sha256').write_text('a' * 64 + ' fixture.sif\n')
            with (validation / 'r/validated_code.tsv').open('w') as stream:
                writer = csv.writer(stream, delimiter='\t')
                writer.writerow(['file', 'md5'])
                for name in ('analysis.R', 'prepare.py', 'model_worker.R', 'recover_checkpoint.R', 'run_dmrcate.sh'):
                    writer.writerow([name, '0' * 32 if name == 'run_dmrcate.sh' else hashlib.md5((WORKFLOW / name).read_bytes()).hexdigest()])
            config = root / 'config.json'
            config.write_text(json.dumps(dict(runtime_validation=str(validation), runtime=str(runtime))))
            result = subprocess.run(['bash', '-c', prefix, 'fixture', str(config), '2026-11-02_13-00-fixture-v1'],
                env=dict(os.environ, DMRCATE_WORKFLOW=str(WORKFLOW)), text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Current code differs from validated fixtures', result.stderr)
            self.assertIn('run_dmrcate.sh', result.stderr)
            self.assertEqual(set(p.name for p in root.iterdir()), {'validation', 'runtime', 'config.json'})


if __name__ == '__main__':
    unittest.main()
