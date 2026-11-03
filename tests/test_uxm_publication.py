"""Exercise the standalone launcher with mocked PAT/UXM tools and tiny files."""
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import shlex
import tempfile
import time
import unittest


class UXMPublicationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir='/tmp')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.repo = self.root / 'repo'
        self.workflow = self.repo / 'workflows/uxm'
        self.workflow.mkdir(parents=True)
        self.results = self.root / 'results'; self.results.mkdir()
        self.workbase = self.root / 'work'
        (self.workflow / 'results').symlink_to(self.results)
        self.bin = self.root / 'bin'; self.bin.mkdir()
        self.executable(self.bin/'python3', '#!/bin/sh\nif [ "$1" = -c ]; then exit 0; fi\nexec '+shlex.quote(sys.executable)+' "$@"\n')
        (self.workflow / 'software/UXM_deconv/.git').mkdir(parents=True)
        atlas = self.workflow / 'reference/Atlas.U25.l4.hg38.full.tsv'
        atlas.parent.mkdir(); atlas.write_text('fixture atlas')
        source = Path(__file__).resolve().parents[1]
        (self.repo / 'scripts/lib').mkdir(parents=True)
        for filename in ('result_names.sh', 'output_guard.py'):
            (self.repo / 'scripts/lib' / filename).write_text((source / 'scripts/lib' / filename).read_text())
        (self.repo/'workflows/output_contract.py').write_text((source/'workflows/output_contract.py').read_text())
        (self.repo / 'metadata').mkdir()
        (self.repo / 'metadata/metadata.tsv').write_text('sample\nS1\n')
        self.input = self.root / 'input'; self.input.mkdir()
        (self.input / 'pat_provenance.json').write_text('{}')
        prepare = self.repo / 'workflows/pat/prepare_pat.py'
        prepare.parent.mkdir(parents=True)
        self.executable(prepare, '''#!/usr/bin/env python3
import sys
from pathlib import Path
p=Path(sys.argv[sys.argv.index('--output-dir')+1]); p.mkdir()
(p/'pat').mkdir()
(p/'pat/S1.pat.gz').write_text('pat')
(p/'pat/S1.pat.gz.csi').write_text('index')
(p/'pat_manifest.tsv').write_text('sample\\tpat\\nS1\\tpat/S1.pat.gz\\n')
for n in ('pat_provenance.json','pat_execution.json','PAT_STARTED','PAT_COMPLETED'):
    (p/n).write_text('{}')
''')
        conda = self.root / 'conda'
        (conda / 'etc/profile.d').mkdir(parents=True)
        (conda / 'etc/profile.d/conda.sh').write_text('conda() { return 0; }\n')
        self.executable(self.bin / 'conda', f'#!/bin/sh\necho "{conda}"\n')
        self.executable(self.bin / 'git', '#!/bin/sh\ncase "$*" in\n  *status*) exit 0;;\n  *--show-toplevel*) echo '+shlex.quote(str(self.repo))+';;\n  *) echo 8d0bb456742ece21802bc8e651e080869c67c944;;\nesac\n')
        self.executable(self.bin / 'sha256sum', '#!/bin/sh\ncat >/dev/null\n')
        for name in ('samtools', 'tabix', 'wgbstools'):
            self.executable(self.bin / name, '#!/bin/sh\necho fixture-tool\n')
        self.executable(self.bin / 'uxm', '''#!/usr/bin/env python3
import os,sys,time
from pathlib import Path
if sys.argv[1]=='deconv':
    Path(os.environ['MOCK_ENTERED']).touch()
    while os.environ.get('MOCK_RELEASE') and not Path(os.environ['MOCK_RELEASE']).exists():
        time.sleep(.02)
    if os.environ.get('MOCK_COLLISION'):
        p=Path(os.environ['MOCK_FINAL'])/'results'; p.mkdir(); (p/'sentinel').write_text('original')
    if os.environ.get('MOCK_FAIL'): sys.exit(7)
    Path(sys.argv[sys.argv.index('--output')+1]).write_text('cell,S1\\nB,1\\n')
elif sys.argv[1]=='plot':
    Path(sys.argv[sys.argv.index('--outpath')+1]).write_text('pdf')
else: print('fixture-uxm')
''')
        launcher = (source / 'workflows/uxm/run_lineage_deconvolution.sh').read_text()
        launcher = launcher.replace('/path/to/hpc-home/analysis/CIRCALL_methylomics', str(self.repo)).replace(
            '/path/to/hpc-scratch/tmp/em_seq/downstream/lineage_deconvolution', str(self.workbase)).replace(
            '/path/to/hpc-scratch/nextflow/nf-core-methylseq/results/uxm_deconvolution', str(self.results))
        self.launcher = self.workflow / 'run_lineage_deconvolution.sh'
        self.launcher.write_text(launcher)
        self.name = '2026-10-28_14-35-uxm-fixture-v1'
        self.final = self.results / (self.name + '-uxm-full-atlas')
        self.work = self.workbase / ('.' + self.final.name + '.running-123')
        self.env = dict(os.environ, SLURM_JOB_ID='123', UXM_RUN_BASE_NAME=self.name,
                        PAT_INPUT_DIR=str(self.input), METHYLSEQ_RESULTS_DIR=str(self.input),
                        PATH=str(self.bin)+':'+os.environ['PATH'], MOCK_FINAL=str(self.final),
                        MOCK_ENTERED=str(self.root/'entered'))

    def executable(self, path, text):
        path.write_text(text); path.chmod(0o755)

    def run_launcher(self, **env):
        return subprocess.run(['bash', str(self.launcher)], env=dict(self.env, **env),
                              capture_output=True, text=True, timeout=20)

    def test_success_publishes_and_records_execution(self):
        result = self.run_launcher()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(self.work.is_dir())
        self.assertEqual(list(self.work.iterdir()), [])
        record = json.loads((self.final/'execution.json').read_text())
        self.assertEqual(record['status'], 'completed')
        self.assertLessEqual(datetime.fromisoformat(record['started_utc'].replace('Z', '+00:00')),
                             datetime.fromisoformat(record['completed_utc']))
        self.assertTrue(Path(str(self.final)+'.lock').exists())

    def test_existing_file_directory_and_dangling_link_refused(self):
        for kind in ('file', 'directory', 'symlink'):
            with self.subTest(kind=kind):
                if kind == 'file': self.final.write_text('original')
                elif kind == 'directory': self.final.mkdir()
                else: self.final.symlink_to(self.root/'absent')
                result = self.run_launcher()
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.work.exists())
                if kind == 'directory': self.final.rmdir()
                else:
                    if kind == 'file': self.assertEqual(self.final.read_text(), 'original')
                    else: self.assertTrue(self.final.is_symlink())
                    self.final.unlink()

    def test_late_collision_retains_work_and_failure_timing(self):
        result = self.run_launcher(MOCK_COLLISION='1')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((self.final/'results/sentinel').read_text(), 'original')
        record = json.loads((self.work/'execution.json').read_text())
        self.assertEqual(record['status'], 'failed_or_interrupted')
        self.assertIn('ended_utc', record)
        self.assertNotIn('completed_utc', record)
        self.assertTrue((self.work/'results/lineage_fractions.csv').exists())

    def test_analysis_failure_records_no_success(self):
        result = self.run_launcher(MOCK_FAIL='1')
        self.assertEqual(result.returncode, 7, result.stderr)
        record = json.loads((self.work/'execution.json').read_text())
        self.assertIn('started_utc', record)
        self.assertNotIn('completed_utc', record)
        self.assertEqual(json.loads((self.final/'.run-ownership.json').read_text())['status'], 'failed')

    def test_concurrent_submission_cannot_reserve_work(self):
        release = self.root/'release'
        first = subprocess.Popen(['bash', str(self.launcher)], env=dict(self.env, MOCK_RELEASE=str(release)),
                                 stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        def stop_first():
            if first.poll() is None: first.kill()
            first.wait()
            first.stderr.close()
        self.addCleanup(stop_first)
        deadline = time.monotonic()+15
        while not (self.root/'entered').exists() and first.poll() is None and time.monotonic()<deadline:
            time.sleep(.02)
        self.assertTrue((self.root/'entered').exists())
        second = self.run_launcher(SLURM_JOB_ID='124')
        self.assertNotEqual(second.returncode, 0)
        self.assertRegex(second.stderr, 'locked|already exists')
        self.assertFalse((self.workbase/('.'+self.final.name+'.running-124')).exists())
        release.touch()
        _, stderr = first.communicate(timeout=15)
        self.assertEqual(first.returncode, 0, stderr)
        self.assertTrue(self.final.is_dir())

    def test_invalid_names_refused_before_work(self):
        for name in ('2026-02-30_14-35-bad-v1', '2026-10-28_24-00-bad-v1', '2026-10-28-bad-v1'):
            with self.subTest(name=name):
                self.assertNotEqual(self.run_launcher(UXM_RUN_BASE_NAME=name).returncode, 0)
                self.assertFalse(self.workbase.exists())


if __name__ == '__main__':
    unittest.main()
