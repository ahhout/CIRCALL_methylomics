"""Mock Slurm/software commands while checking dispatch and the shared lock."""
from pathlib import Path
import os
import shlex
import subprocess
import tempfile
import unittest

class ShellDispatch(unittest.TestCase):
    def test_each_method_dispatches_once_and_holds_a_fixture_lock(self):
        source = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory(dir='/tmp') as directory:
            root = Path(directory); tools = root/'bin'; tools.mkdir()
            conda = root/'conda'; (conda/'etc/profile.d').mkdir(parents=True)
            (conda/'etc/profile.d/conda.sh').write_text('conda() { :; }; module() { :; }\n')
            def executable(name, body):
                p=tools/name; p.write_text('#!/bin/bash\n'+body); p.chmod(0o755)
            executable('conda', 'echo '+shlex.quote(str(conda))+'\n')
            executable('apptainer', 'echo "apptainer version 1.3.6"\n')
            for name in ('nextflow','sbatch'): executable(name, 'exit 0\n')
            executable('python', 'test -e /proc/$$/fd/8 || exit 9\nprintf "%s\\n" "$*" >> "$DISPATCH_LOG"\n')
            launcher=root/'run.sh'
            launcher.write_text((source/'run_deconflow.sh').read_text().replace(
                '/path/to/hpc-scratch/tmp/em_seq/.methylseq-cache-maintenance.lock', str(root/'cache.lock')))
            env=dict(os.environ, PATH=str(tools)+':'+os.environ['PATH'], SLURM_JOB_ID='123',
                     DECONFLOW_WORKFLOW=str(source),DISPATCH_LOG=str(root/'dispatch.txt'))
            for method, expected in (('uxm','analysis.py'),('celfie','tools/celfie_patients.py')):
                result=subprocess.run(['bash',str(launcher)]+(['celfie'] if method=='celfie' else [])+['fixture-output'],
                                      env=env,cwd='/tmp',capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                lines=(root/'dispatch.txt').read_text().splitlines()
                self.assertEqual(lines[-1],str(source/expected)+' fixture-output')
            self.assertEqual(len(lines),2)

if __name__ == '__main__': unittest.main()
