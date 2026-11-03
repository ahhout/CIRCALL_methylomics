import copy
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('cache',Path(__file__).resolve().parents[1]/'scripts/cleanup_methylseq_cache.py')
cache=importlib.util.module_from_spec(spec);spec.loader.exec_module(cache)


class CacheTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir='/tmp');self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name)
        self.roots=tuple(self.base/n for n in ('integrated','old','spike','atlas'))
        for r in self.roots:r.mkdir()
        self.patch=patch.multiple(cache,ROOTS=self.roots,BASE=self.base,LOCK=self.base/'lock',
                                  ATLAS_LOCK=self.base/'atlas.lock')
        self.patch.start();self.addCleanup(self.patch.stop)
        self.task=self.make_task(self.roots[1],'aa','1'*30)

    def make_task(self,root,prefix,digest):
        t=root/prefix/digest;t.mkdir(parents=True)
        (t/'payload.bam').write_bytes(b'payload')
        (t/'.command.sh').write_text('echo task\n')
        return t

    def state(self):
        for root in self.roots:
            for p in root.rglob('*'):
                if p.is_file() and not p.is_symlink():
                    with p.open('rb') as h: os.fsync(h.fileno())
        scans=[cache.inspect(r,True) for r in self.roots]
        tasks,candidates,protected,errors,findings=cache.assess(scans,[],{})
        return dict(tasks=tasks,candidates=candidates,protected=sorted(protected),errors=errors,
                    findings=findings,traces=[],queue={'jobs':[],'error':None})

    def test_protected_dependency_closure(self):
        e31=self.make_task(self.roots[0],'bb','2'*30)
        (e31/'input.bam').symlink_to(self.task/'payload.bam')
        scans=[cache.inspect(r,True) for r in self.roots]
        _,candidates,protected,_,_=cache.assess(scans,[],{str(e31):set(cache.PROTECTED_RUNS)})
        self.assertIn(str(self.task),protected);self.assertNotIn(str(self.task),candidates)

    def test_container_mount_is_not_a_dependency_on_every_task(self):
        e31=self.make_task(self.roots[0],'bb','2'*30)
        script=e31/'.command.run'
        script.write_text('apptainer exec -B '+str(self.roots[1])+' image.sif bash script\n')
        scans=[cache.inspect(r,True) for r in self.roots]
        _,cand,_,_,findings=cache.assess(scans,[],{str(e31):set(cache.PROTECTED_RUNS)})
        self.assertIn(str(self.task),cand)
        self.assertTrue(any(f['kind']=='container_mount' for f in findings))
        with script.open('a') as h:h.write('cat '+str(self.task/'payload.bam')+'\n')
        _,cand,protected,_,_=cache.assess([cache.inspect(r,True) for r in self.roots],[],
                                         {str(e31):set(cache.PROTECTED_RUNS)})
        self.assertIn(str(self.task),protected)

    def test_non_mount_root_reference_keeps_whole_root_protected(self):
        e31=self.make_task(self.roots[0],'bb','2'*30)
        (e31/'.command.run').write_text('apptainer exec -B '+str(self.roots[1])+\
            ' image.sif bash script\nfind '+str(self.roots[1])+'\n')
        _,_,protected,_,_=cache.assess([cache.inspect(r,True) for r in self.roots],[],
                                      {str(e31):set(cache.PROTECTED_RUNS)})
        self.assertIn(str(self.task),protected)

    def test_chained_external_link_protects_target(self):
        ext=self.base/'external';ext.mkdir()
        (ext/'one').symlink_to(self.task/'payload.bam');(ext/'two').symlink_to(ext/'one')
        _,cand,protected,_,_=cache.assess([cache.inspect(r,True) for r in self.roots],[cache.inspect(ext)],{})
        self.assertIn(str(self.task),protected)

    def test_nextflow_config_reference_requires_review(self):
        ext=self.base/'external';ext.mkdir()
        (ext/'nextflow.config').write_text('params.fasta = "'+str(self.task/'payload.bam')+'"\n')
        _,_,_,_,findings=cache.assess([cache.inspect(r,True) for r in self.roots],
                                     [cache.inspect(ext)],{})
        self.assertTrue(any(f['kind']=='text_reference' and str(self.task) in f['tasks'] for f in findings))

    def test_intermediate_task_in_symlink_chain_is_protected(self):
        middle=self.make_task(self.roots[3],'ee','5'*30)
        (middle/'via').symlink_to(self.task/'payload.bam')
        protected=self.make_task(self.roots[0],'bb','2'*30)
        (protected/'input').symlink_to(middle/'via')
        _,_,kept,_,_=cache.assess([cache.inspect(r,True) for r in self.roots],[],
                                 {str(protected):set(cache.PROTECTED_RUNS)})
        self.assertIn(str(middle),kept);self.assertIn(str(self.task),kept)

    def test_broken_link_blocks_review(self):
        (self.task/'bad').symlink_to(self.roots[0]/'missing-input')
        self.assertTrue(any('Unresolved link' in e for e in self.state()['errors']))

    def test_only_reviewed_missing_external_link_can_be_resolved(self):
        ext=self.base/'external';ext.mkdir()
        link=ext/'historical.pat.gz'
        link.symlink_to('/path/to/hpc-home/nonexistent-retired-pilot/path.pat.gz')
        _,_,_,errors,findings=cache.assess([cache.inspect(r,True) for r in self.roots],
                                         [cache.inspect(ext)],{})
        state={'errors':errors,'findings':findings}
        self.assertTrue(cache.unresolved_errors(state,None))
        self.assertEqual(findings[0]['kind'],'missing_external_link')
        review={'findings':{cache.sha(cache.encode(findings[0])):
                 {'status':'historical','reason':'fixture historical provenance inspected'}}}
        self.assertFalse(cache.unresolved_errors(state,review))
        (ext/'cache-input').symlink_to(self.roots[0]/'missing-input')
        _,_,_,errors,findings=cache.assess([cache.inspect(r,True) for r in self.roots],
                                         [cache.inspect(ext)],{})
        self.assertTrue(cache.unresolved_errors({'errors':errors,'findings':findings},review))

    def test_broken_external_chain_through_cache_cannot_be_waived(self):
        ext=self.base/'external';ext.mkdir()
        (self.task/'missing').symlink_to('/path/to/hpc-home/nonexistent-retired-pilot/path')
        (ext/'consumer').symlink_to(self.task/'missing')
        _,candidates,protected,errors,findings=cache.assess(
            [cache.inspect(r,True) for r in self.roots],[cache.inspect(ext)],{})
        self.assertIn(str(self.task),protected)
        self.assertNotIn(str(self.task),candidates)
        self.assertFalse(any(f['kind']=='missing_external_link' for f in findings))
        self.assertTrue(errors)

    def test_metadata_scan_has_bounded_longer_timeout(self):
        import subprocess
        with patch.object(cache.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'{"_scan_errors": []}\n')) as run:
            cache.scan(self.roots[0],True)
            self.assertEqual(run.call_args.kwargs['timeout'],900)
            cache.scan(self.base,False)
            self.assertEqual(run.call_args.kwargs['timeout'],900)

    def test_hardlink_excludes_task(self):
        import os
        os.link(self.task/'payload.bam',self.base/'outside.bam')
        self.assertNotIn(str(self.task),self.state()['candidates'])

    def test_unmatched_integrated_task_is_retained(self):
        new=self.make_task(self.roots[0],'dd','4'*30)
        self.assertIn(str(new),self.state()['protected'])

    def test_other_batch_dependency_prevents_retirement(self):
        other=self.make_task(self.roots[3],'cc','3'*30)
        (other/'input.bam').symlink_to(self.task/'payload.bam')
        self.assertNotIn(str(self.task),cache.batch_selection(self.state(),'legacy_emseq'))

    def test_active_and_unavailable_scheduler_refused(self):
        for queue in ({'error':'unavailable','jobs':[]},{'error':None,'jobs':['123|nf_methylseq|/work|cmd']}):
            with self.assertRaises(cache.Refused):cache.require_idle({},queue)

    def test_independent_atlas_only_accepts_exact_reviewed_unrelated_jobs(self):
        job='42|nf_methylseq|'+str(self.roots[0])+'|.command.run'
        review={'maintenance_policy':'independent_atlas',
                'concurrent_job_evidence':{cache.sha(job.encode()):'Reviewed integrated inputs'}}
        cache.require_idle(review,{'error':None,'jobs':[job]})
        for bad in (job.replace('42','43',1),'44|bcp_all_wgbs|/work|script'):
            review['concurrent_job_evidence'][cache.sha(bad.encode())]='Cannot waive atlas jobs'
            if 'bcp_all' in bad:
                with self.assertRaises(cache.Refused):cache.require_idle(review,{'error':None,'jobs':[bad]})
            else:
                del review['concurrent_job_evidence'][cache.sha(bad.encode())]
                with self.assertRaises(cache.Refused):cache.require_idle(review,{'error':None,'jobs':[bad]})
        with self.assertRaises(cache.Refused):cache.require_idle({}, {'error':None,'jobs':[job]})

    def test_atlas_lock_blocks_cleanup_but_not_integrated_lock(self):
        with cache.ATLAS_LOCK.open('a') as h:
            fcntl.flock(h,fcntl.LOCK_SH|fcntl.LOCK_NB)
            with self.assertRaises(BlockingIOError):
                with cache.maintenance(cache.ATLAS_LOCK):pass
            with cache.maintenance():pass

    def test_independent_review_keeps_atlas_and_unscanned_findings(self):
        findings=[{'path':'/log','kind':'text_reference','target':str(self.roots[0])},
                  {'path':'/log','kind':'text_reference','target':str(self.roots[3])},
                  {'path':'/large','kind':'oversized_text'}]
        self.assertEqual(cache.review_findings({'findings':findings},
                         {'maintenance_policy':'independent_atlas'}),findings[1:])

    def test_independent_policy_refuses_other_batch(self):
        path,state=self.manifest();m=json.loads(path.read_text())
        m['review']['maintenance_policy']='independent_atlas';path.write_bytes(cache.encode(m))
        with self.assertRaisesRegex(cache.Refused,'restricted'):
            cache.apply(path,cache.sha(path.read_bytes()))

    def test_independent_atlas_apply_allows_unrelated_output_growth(self):
        path,_=self.manifest();m=json.loads(path.read_text())
        atlas=self.make_task(self.roots[3],'ee','5'*30)
        active=self.make_task(self.roots[0],'bb','2'*30)
        before=self.state();selected=[str(atlas)]
        run=sorted(cache.PROTECTED_RUNS)[-1]
        before['memberships']={str(active):[run]}
        job='42|nf_methylseq|'+str(active)+'|.command.run'
        m.update(batch='atlas_wgbs',state=before,tasks=selected,entries=before['tasks'][str(atlas)])
        m['review'].update(maintenance_policy='independent_atlas',concurrent_runs=[run],tasks={str(atlas):
            {'no_active_consumers':True,'evidence':'Independent atlas fixture'}},
            atlas_published_inputs_independent=True,concurrent_job_evidence={cache.sha(job.encode()):'Reviewed fixture'})
        path.write_bytes(cache.encode(m))
        (active/'payload.bam').write_bytes(b'active data grew')
        current=self.state();post=copy.deepcopy(current)
        del post['tasks'][str(atlas)];del post['candidates'][str(atlas)]
        with patch.object(cache,'collect',side_effect=[current,post]), \
             patch.object(cache,'scheduler',return_value={'jobs':[job],'error':None}), \
             patch.object(cache,'scan',side_effect=lambda p,c=False:cache.inspect(p,c)):
            cache.apply(path,cache.sha(path.read_bytes()))
        self.assertFalse(atlas.exists());self.assertEqual((active/'payload.bam').read_bytes(),b'active data grew')
        self.assertTrue((path.parent/'execution.json').exists())

    def test_independent_atlas_still_detects_missing_protected_directory(self):
        active=self.make_task(self.roots[0],'bb','2'*30);before=self.state();after=copy.deepcopy(before)
        del after['tasks'][str(active)]
        with self.assertRaises(cache.Refused):cache.check_protected(before,after,{'maintenance_policy':'independent_atlas'})

    def test_independent_atlas_keeps_canonical_file_checks(self):
        e31=self.make_task(self.roots[0],'bb','2'*30);before=self.state()
        before['memberships']={str(e31):[sorted(cache.PROTECTED_RUNS)[0]]}
        (e31/'payload.bam').write_bytes(b'unexpected canonical change')
        review={'maintenance_policy':'independent_atlas','concurrent_runs':[sorted(cache.PROTECTED_RUNS)[-1]]}
        with self.assertRaises(cache.Refused):cache.check_protected(before,self.state(),review)

    def test_independent_atlas_excludes_live_input_dependencies(self):
        atlas=self.make_task(self.roots[3],'ee','5'*30)
        active=self.make_task(self.roots[0],'bb','2'*30)
        (active/'input').symlink_to(atlas/'payload.bam')
        self.assertNotIn(str(atlas),cache.batch_selection(self.state(),'atlas_wgbs'))

    def test_maintenance_lock_contention(self):
        with cache.maintenance():
            with self.assertRaises(BlockingIOError):
                with cache.maintenance():pass

    def test_ambiguous_trace_hash(self):
        other=self.make_task(self.roots[1],'aa','111111'+'2'*24)
        results=self.base/'results';d=results/'bulk'/next(iter(cache.PROTECTED_RUNS))/'pipeline_info';d.mkdir(parents=True)
        (d/'execution_trace.txt').write_text('hash\tname\naa/111111\ttask\n')
        with patch.object(cache,'RESULTS',results):
            _,errors,_=cache.memberships([str(self.task),str(other)])
        self.assertTrue(any('Ambiguous' in e for e in errors))

    def test_interrupted_attempt_membership_from_full_log_path(self):
        results=self.base/'results';run=next(iter(cache.PROTECTED_RUNS))
        d=results/'bulk'/run/'pipeline_info';d.mkdir(parents=True)
        (d/'nextflow_123.log').write_text('jobId: 123; workDir: '+str(self.task)+'\n')
        with patch.object(cache,'RESULTS',results):
            members,_,evidence=cache.memberships([str(self.task)])
        self.assertIn(run,members[str(self.task)])
        self.assertEqual(evidence[0]['kind'],'saved_task_log')

    def test_atlas_final_and_sensitivity_traces_are_included(self):
        results=self.base/'results'
        for mode in ('final','sensitivity'):
            d=results/'bcp_all_atlas_step3'/mode/('atlas_'+mode)/'pipeline_info'
            d.mkdir(parents=True)
            (d/'execution_trace.txt').write_text('hash\tname\naa/111111\tatlas task\n')
        with patch.object(cache,'RESULTS',results):
            members,_,_=cache.memberships([str(self.task)])
        self.assertEqual(members[str(self.task)],{'atlas_final','atlas_sensitivity'})

    def test_changed_or_added_entry_refused(self):
        entries=self.state()['tasks'][str(self.task)]
        (self.task/'new').write_text('new')
        with patch.object(cache,'scan',side_effect=lambda p,c=False:cache.inspect(p,c)):
            with self.assertRaisesRegex(cache.Refused,'changed'):cache.verify_entries(entries,[str(self.task)])

    def test_symlink_target_is_never_unlinked(self):
        target=self.base/'outside';target.write_text('keep')
        link=self.task/'input';link.symlink_to(target)
        row=dict(path=str(link),**cache.identity(link));cache.remove_entry(row)
        self.assertEqual(target.read_text(),'keep');self.assertFalse(link.is_symlink())

    def test_replaced_directory_symlink_refused(self):
        target=self.base/'outside';target.mkdir()
        directory=self.task/'sub';directory.mkdir();payload=directory/'file';payload.write_text('old')
        row=dict(path=str(payload),**cache.identity(payload))
        payload.unlink();directory.rmdir();directory.symlink_to(target,target_is_directory=True)
        (target/'file').write_text('keep')
        with self.assertRaises(OSError):cache.remove_entry(row)
        self.assertEqual((target/'file').read_text(),'keep')

    def manifest(self):
        state=self.state();selected=cache.batch_selection(state,'legacy_emseq')
        evidence=self.base/'evidence.json';evidence.write_text('{"retained": true}')
        review=dict(checked_utc=cache.now(),tasks={t:{'no_active_consumers':True,'evidence':'fixture'} for t in selected},
                    recovery_policy='retire_resume',findings={},retained_evidence=[{'path':str(evidence),'sha256':cache.sha(evidence.read_bytes())}])
        m=dict(tasks=selected,entries=[r for t in selected for r in state['tasks'][t]],state=state,
               batch='legacy_emseq',review=review,blockers=[],helper_sha256=cache.sha(Path(cache.__file__).read_bytes()))
        out=self.base/'audit';out.mkdir();manifest=out/'manifest.json';manifest.write_bytes(cache.encode(m))
        return manifest,state

    def test_manifest_edit_after_approval_is_refused(self):
        path,state=self.manifest();approved=cache.sha(path.read_bytes())
        path.write_bytes(path.read_bytes()+b'\n')
        with self.assertRaisesRegex(cache.Refused,'hash mismatch'):cache.apply(path,approved)

    def test_stale_review_is_refused(self):
        path,state=self.manifest();m=json.loads(path.read_text())
        m['review']['checked_utc']='2000-01-01T00:00:00+00:00'
        self.assertTrue(cache.validate_review(m['review'],state,m['tasks']))

    def test_added_file_after_approval_blocks_apply(self):
        path,state=self.manifest();(self.task/'unexpected').write_text('keep')
        with patch.object(cache,'collect',return_value=state), \
             patch.object(cache,'scheduler',return_value={'jobs':[],'error':None}), \
             patch.object(cache,'scan',side_effect=lambda p,c=False:cache.inspect(p,c)):
            with self.assertRaisesRegex(cache.Refused,'changed'):cache.apply(path,cache.sha(path.read_bytes()))
        self.assertTrue((self.task/'payload.bam').exists())

    def test_wrong_approval_hash_refused(self):
        path,state=self.manifest()
        with self.assertRaisesRegex(cache.Refused,'hash mismatch'):cache.apply(path,'0'*64)
        self.assertTrue((self.task/'payload.bam').exists())

    def test_successful_apply_preserves_provenance_and_journals(self):
        path,state=self.manifest()
        with patch.object(cache,'collect',side_effect=[state,dict(state,tasks={},candidates={})]), \
             patch.object(cache,'scheduler',return_value={'jobs':[],'error':None}), \
             patch.object(cache,'scan',side_effect=lambda p,c=False:cache.inspect(p,c)):
            cache.apply(path,cache.sha(path.read_bytes()))
        self.assertFalse(self.task.exists())
        journal=[json.loads(l) for l in (path.parent/'deletions.jsonl').read_text().splitlines()]
        self.assertEqual(sum(r['event']=='intent' for r in journal),sum(r['event']=='removed' for r in journal))
        self.assertTrue(list((path.parent/'task_provenance').rglob('.command.sh')))

    def test_interrupted_apply_reports_partial_completion(self):
        path,state=self.manifest();original=cache.remove_entry;calls=[]
        def fail_second(row):
            calls.append(row)
            if len(calls)==2:raise cache.Refused('injected stop')
            original(row)
        with patch.object(cache,'collect',return_value=state), \
             patch.object(cache,'scheduler',return_value={'jobs':[],'error':None}), \
             patch.object(cache,'scan',side_effect=lambda p,c=False:cache.inspect(p,c)), \
             patch.object(cache,'remove_entry',side_effect=fail_second):
            with self.assertRaisesRegex(cache.Refused,'injected'):cache.apply(path,cache.sha(path.read_bytes()))
        events=[json.loads(l)['event'] for l in (path.parent/'deletions.jsonl').read_text().splitlines()]
        self.assertEqual(events[-1],'stopped');self.assertNotIn('complete',events)
        self.assertFalse((path.parent/'execution.json').exists())


if __name__=='__main__':unittest.main()
