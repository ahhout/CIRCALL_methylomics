import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

# Frozen historical entrypoint is tested only inside temporary fake runtime trees.
LEGACY_SOURCE = Path(__file__).resolve().parents[1] / 'archive/preprocessing/legacy_project_entrypoint'
spec = importlib.util.spec_from_file_location('session', LEGACY_SOURCE / 'scripts/methylseq_cache_session.py')
session = importlib.util.module_from_spec(spec); spec.loader.exec_module(session)
SID = '17da6065-5303-40ef-8047-5c8ff69a052d'
OTHER = 'aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee'


class SessionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir='/tmp'); self.addCleanup(self.tmp.cleanup)
        self.repo = Path(self.tmp.name)
        (self.repo / '.nextflow/cache' / SID / 'db').mkdir(parents=True)
        (self.repo / '.nextflow/cache' / SID / 'db/CURRENT').write_text('MANIFEST-test\n')
        self.record = self.repo / 'record.json'
        self.value = dict(schema=1, session_id=SID, work_dir=session.EXPECTED_WORK)
        self.record.write_text(json.dumps(self.value))
        self.history = self.repo / '.nextflow/history'
        self.history.write_text(f'date\tduration\tname\tOK\tcommit\t{SID}\tnextflow run pipeline -work-dir {session.EXPECTED_WORK} --outdir /example\n')

    def test_unrelated_latest_history_does_not_change_selection(self):
        with self.history.open('a') as f:
            f.write(f'date\t-\tother\tOK\tx\t{OTHER}\tnextflow run other -work-dir /other\n')
        self.assertEqual(session.validate_session(self.record, self.repo)['session_id'], SID)

    def test_missing_cache_fails(self):
        (self.repo / '.nextflow/cache' / SID / 'db/CURRENT').unlink()
        with self.assertRaisesRegex(ValueError, 'missing'):
            session.validate_session(self.record, self.repo)

    def test_conflicting_history_fails(self):
        with self.history.open('a') as f:
            f.write(f'd\t-\tn\tOK\tx\t{OTHER}\tnextflow run p -work-dir {session.EXPECTED_WORK}\n')
        with self.assertRaisesRegex(ValueError, 'conflicting'):
            session.validate_session(self.record, self.repo)

    def test_missing_record_fails(self):
        self.record.unlink()
        with self.assertRaises(OSError): session.validate_session(self.record, self.repo)

    def snapshots(self):
        submitted = self.repo / 'submitted'; submitted.mkdir()
        paths = {'params.emseq.integrated.yml': 'config/params.emseq.integrated.yml',
                 'nextflow.config': 'config/nextflow.config', 'samplesheet.csv': 'samplesheets/samplesheet.csv',
                 'reference_manifest.tsv': 'references/emseq_hg38_lambda_puc19/manifest.tsv'}
        names = dict(paths, **{'run_methylseq_integrated.sh':'run_methylseq_integrated.sh'})
        sums=[]
        for name, relative in names.items():
            data=('original '+name).encode()
            (submitted/name).write_bytes(data)
            path=self.repo/relative; path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
            sums.append(hashlib.sha256(data).hexdigest()+'  '+name)
        (submitted/'SHA256SUMS').write_text('\n'.join(sums)+'\n')
        return submitted

    def test_launcher_revision_does_not_change_scientific_input_contract(self):
        submitted=self.snapshots()
        (self.repo/'run_methylseq_integrated.sh').write_text('new bookkeeping')
        session.verify_snapshots(submitted,self.repo)
        self.assertEqual(session.infer_output_session(Path('/example'),self.repo),SID)

    def test_changed_scientific_input_fails(self):
        submitted=self.snapshots()
        (self.repo/'samplesheets/samplesheet.csv').write_text('changed samples')
        with self.assertRaisesRegex(ValueError,'scientific input'):
            session.verify_snapshots(submitted,self.repo)

    def test_tampered_original_launcher_fails(self):
        submitted=self.snapshots();(submitted/'run_methylseq_integrated.sh').write_text('tampered')
        with self.assertRaisesRegex(ValueError,'snapshot changed'):
            session.verify_snapshots(submitted,self.repo)

    def test_ambiguous_legacy_output_mapping_fails(self):
        with self.history.open('a') as h:
            h.write(f'd\t-\tn\tOK\tx\t{OTHER}\tnextflow run p --outdir /example\n')
        with self.assertRaisesRegex(ValueError,'unambiguous'):
            session.infer_output_session(Path('/example'),self.repo)

    def test_launcher_uses_explicit_uuid_and_exit_recording(self):
        text=(LEGACY_SOURCE/'run_methylseq_integrated.sh').read_text()
        self.assertIn('-resume "$SESSION_ID"',text)
        self.assertNotIn('    -resume\n',text)
        self.assertIn('trap record_attempt_exit EXIT',text)
        self.assertIn('flock -n 9',text)

    def test_interrupted_attempt_records_incomplete_membership(self):
        output=self.repo/'out';attempt=output/'pipeline_info/submitted/attempts/1';attempt.mkdir(parents=True)
        (attempt/'provenance.json').write_text('{}')
        session.finish(output,'1',143)
        record=json.loads((attempt/'provenance.json').read_text())
        self.assertFalse(record['membership_complete'])
        self.assertEqual(record['status'],'failed_or_interrupted')



class LauncherIntegrationTests(unittest.TestCase):
    """Execute the launcher with mocked cluster/runtime commands, never methylseq."""
    def setUp(self):
        import os
        self.tmp=tempfile.TemporaryDirectory(dir='/tmp');self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.repo=self.root/'repo';self.repo.mkdir()
        self.work=self.root/'scratch/work_emseq_integrated';self.work.mkdir(parents=True)
        self.results=self.root/'results';self.bin=self.root/'bin';self.bin.mkdir()
        source=LEGACY_SOURCE
        for sub in ('scripts','config','samplesheets','references/emseq_hg38_lambda_puc19/BismarkIndex/Bisulfite_Genome/CT_conversion',
                    'references/emseq_hg38_lambda_puc19/BismarkIndex/Bisulfite_Genome/GA_conversion','.nextflow/cache/'+SID+'/db'):
            (self.repo/sub).mkdir(parents=True,exist_ok=True)
        (self.repo/'.nextflow/cache'/SID/'db/CURRENT').write_text('MANIFEST-test\n')
        launcher=(source/'run_methylseq_integrated.sh').read_text().replace('/path/to/hpc-home/em_seq',str(self.repo)).replace(
            '/path/to/hpc-scratch/tmp/em_seq',str(self.root/'scratch')).replace(
            '/path/to/hpc-scratch/nextflow/nf-core-methylseq/results',str(self.results))
        self.launcher=self.repo/'run_methylseq_integrated.sh';self.launcher.write_text(launcher)
        helper=(source/'scripts/methylseq_cache_session.py').read_text().replace(session.EXPECTED_WORK,str(self.work))
        (self.repo/'scripts/methylseq_cache_session.py').write_text(helper)
        for name in ('params.emseq.integrated.yml','nextflow.config'):(self.repo/'config'/name).write_text('original\n')
        (self.repo/'config/methylseq_cache_session.json').write_text(json.dumps(dict(schema=1,session_id=SID,work_dir=str(self.work))))
        (self.repo/'.nextflow/history').write_text(f'd\t-\told\tOK\tx\t{SID}\tnextflow run p -work-dir {self.work}\n'+
            f'd\t-\tunrelated\tOK\tx\t{OTHER}\tnextflow run p -work-dir /unrelated\n')
        ref=self.repo/'references/emseq_hg38_lambda_puc19';(ref/'manifest.tsv').write_text('reference\n')
        (ref/'BismarkIndex/genome.fa').write_text('>phage_lambda\nA\n>plasmid_puc19c\nC\n')
        (ref/'BismarkIndex/Bisulfite_Genome/CT_conversion/index.bt2').write_text('index')
        raw=self.root/'raw.fq';raw.write_text('raw')
        (self.repo/'samplesheets/samplesheet.csv').write_text('sample,fastq_1,fastq_2,genome\nsample,'+str(raw)+',,\n')
        self.executable(self.repo/'scripts/summarize_emseq_controls.py','#!/usr/bin/env python3\n')
        conda=self.root/'conda';(conda/'etc/profile.d').mkdir(parents=True)
        (conda/'etc/profile.d/conda.sh').write_text('conda() { return 0; }\n')
        self.executable(self.bin/'conda',f'#!/bin/sh\necho "{conda}"\n')
        self.executable(self.bin/'apptainer','#!/bin/sh\nexit 0\n')
        self.executable(self.bin/'nextflow',f'''#!/usr/bin/env python3
import sys,os,json
from pathlib import Path
args=sys.argv[1:]
with open({str(self.root/'commands.jsonl')!r},'a') as h:h.write(json.dumps(args)+'\\n')
out=Path(args[args.index('--outdir')+1])
(out/'pipeline_info/execution_trace_mock.txt').write_text('hash\\tname\\tstatus\\n00/abcdef\\tmock\\tCOMPLETED\\n')
sys.exit(int(os.environ.get('MOCK_NF_EXIT','0')))
''')
        envfile=self.root/'bashenv';envfile.write_text('module() { return 0; }\n')
        self.env=dict(os.environ,PATH=str(self.bin)+':'+os.environ['PATH'],BASH_ENV=str(envfile),
                      RUN_NAME='2026-09-21_10-30-emseq-integrated-fixture-v1',SLURM_JOB_ID='1',RESUME_EXISTING='false')
        self.output=self.results/'bulk'/self.env['RUN_NAME']

    def executable(self,p,text):p.write_text(text);p.chmod(0o755)

    def run_launcher(self,**changes):
        import subprocess
        return subprocess.run(['bash',str(self.launcher)],env=dict(self.env,**changes),capture_output=True,text=True,timeout=30)

    def test_new_and_resumed_analysis_pin_same_session(self):
        first=self.run_launcher();self.assertEqual(first.returncode,0,first.stderr)
        self.launcher.write_text(self.launcher.read_text()+'\n# bookkeeping revision\n')
        second=self.run_launcher(SLURM_JOB_ID='2',RESUME_EXISTING='true');self.assertEqual(second.returncode,0,second.stderr)
        commands=[json.loads(l) for l in (self.root/'commands.jsonl').read_text().splitlines()]
        self.assertEqual(len(commands),2)
        for args in commands:self.assertEqual(args[args.index('-resume')+1],SID)
        record=json.loads((self.output/'pipeline_info/submitted/attempts/2/provenance.json').read_text())
        self.assertTrue(record['membership_complete'])
        self.assertIn('started_utc', record)
        self.assertIn('completed_utc', record)
        original=json.loads((self.output/'pipeline_info/submitted/attempts/1/provenance.json').read_text())
        self.assertNotEqual(record['job_id'], original['job_id'])
        self.assertIn('started_utc', original)

    def test_new_outputs_require_valid_hours_and_minutes(self):
        for name in ('2026-09-21-emseq-integrated-fixture-v1',
                     '2026-09-21_24-00-emseq-integrated-fixture-v1',
                     '2026-09-21_10-60-emseq-integrated-fixture-v1'):
            with self.subTest(name=name):
                result = self.run_launcher(RUN_NAME=name)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse((self.root/'commands.jsonl').exists())
                self.assertFalse((self.results/'bulk'/name).exists())

    def test_changed_scientific_inputs_prevent_nextflow_invocation(self):
        self.assertEqual(self.run_launcher().returncode,0)
        (self.repo/'config/params.emseq.integrated.yml').write_text('changed interpretation\n')
        r=self.run_launcher(SLURM_JOB_ID='2',RESUME_EXISTING='true')
        self.assertNotEqual(r.returncode,0)
        self.assertIn('scientific input',r.stderr)
        self.assertEqual(len((self.root/'commands.jsonl').read_text().splitlines()),1)

    def test_failed_attempt_is_marked_incomplete(self):
        r=self.run_launcher(MOCK_NF_EXIT='7');self.assertEqual(r.returncode,7,r.stderr)
        record=json.loads((self.output/'pipeline_info/submitted/attempts/1/provenance.json').read_text())
        self.assertFalse(record['membership_complete']);self.assertEqual(record['exit_status'],7)
        self.assertIn('started_utc', record);self.assertIn('ended_utc', record)
        self.assertNotIn('completed_utc', record)

    def test_maintenance_lock_prevents_nextflow_invocation(self):
        import fcntl
        with (self.root/'scratch/.methylseq-cache-maintenance.lock').open('w') as h:
            fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
            r=self.run_launcher()
        self.assertNotEqual(r.returncode,0);self.assertIn('in use',r.stderr)
        self.assertFalse((self.root/'commands.jsonl').exists())
        self.assertFalse(self.output.exists())

if __name__ == '__main__': unittest.main()
