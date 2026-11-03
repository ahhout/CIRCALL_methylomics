#!/usr/bin/env python3
"""Compare completed CELFIE estimates with the preserved full-atlas UXM baseline."""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from common import ATLAS, ATLAS_SHA, WORKFLOW, atlas_rows, finish, identity, new_run, read_json, rows, sha, utc, write_json
from native import checked


def compare(args):
    _, cells = atlas_rows()
    metadata = rows(WORKFLOW/'config/samples.tsv')
    samples = [r['sample'] for r in metadata]
    source = read_json(WORKFLOW/'config/active_comparison.json')['uxm']['fraction_table']
    uxm_path = Path(source['path'])
    if sha(uxm_path) != source['sha256']:
        raise ValueError('Preserved UXM baseline changed')
    celfie_root = args.celfie.resolve(strict=True)
    gate = celfie_root/'VALIDATED.json'
    if not gate.is_file():
        raise ValueError('Patient CELFIE validation marker is required')
    celfie_path = args.fractions.resolve(strict=True)
    accepted = read_json(gate)
    if accepted['fractions']['sha256'] != sha(celfie_path):
        raise ValueError('CELFIE fractions differ from the validated result')
    uxm = checked(pd.read_csv(uxm_path,index_col=0),samples,cells)
    celfie = checked(pd.read_csv(celfie_path,sep='\t',index_col=0),samples,cells)
    root = new_run(args.output)
    try:
        write_json(root/'inputs.json', {'started_utc':utc(), 'atlas':identity(ATLAS,True),
            'uxm':identity(uxm_path,True), 'celfie':identity(celfie_path,True),
            'validation':identity(gate,True), 'helper':identity(__file__,True)})
        uxm.to_csv(root/'uxm_fractions.tsv',sep='\t',float_format='%.12g')
        celfie.to_csv(root/'celfie_fractions.tsv',sep='\t',float_format='%.12g')
        diff = 100*(celfie-uxm)
        diff.to_csv(root/'celfie_minus_uxm_pp.tsv',sep='\t',float_format='%.12g')
        sample_summary = pd.DataFrame({
            'sum_absolute_difference_pp':diff.abs().sum(axis=1),
            'mean_absolute_difference_pp':diff.abs().mean(axis=1),
            'maximum_absolute_difference_pp':diff.abs().max(axis=1),
            'largest_difference_cell_type':diff.abs().idxmax(axis=1),
            'uxm_largest_cell_type':uxm.idxmax(axis=1),
            'celfie_largest_cell_type':celfie.idxmax(axis=1)})
        sample_summary.index.name='sample'
        sample_summary = pd.DataFrame(metadata).set_index('sample')[['patient','compartment']].join(sample_summary)
        sample_summary.to_csv(root/'sample_summary.tsv',sep='\t',float_format='%.12g')
        cell_summary = pd.DataFrame({
            'mean_uxm_fraction':uxm.mean(axis=0), 'mean_celfie_fraction':celfie.mean(axis=0),
            'mean_difference_pp':diff.mean(axis=0),
            'median_difference_pp':diff.median(axis=0),
            'mean_absolute_difference_pp':diff.abs().mean(axis=0),
            'maximum_absolute_difference_pp':diff.abs().max(axis=0)})
        cell_summary.index.name='cell_type'
        cell_summary.to_csv(root/'cell_type_summary.tsv',sep='\t',float_format='%.12g')
        if not args.no_plots:
            os.environ.setdefault('MPLCONFIGDIR',str(root/'.plot-cache'))
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            colors = plt.cm.gist_ncar(np.linspace(0,.95,len(cells)))
            fig,axes = plt.subplots(2,1,figsize=(16,10),sharex=True)
            for ax,(label,frame) in zip(axes,[('UXM',uxm),('CELFIE: reconstructed megakaryocyte reference',celfie)]):
                frame.mul(100).plot.bar(stacked=True,ax=ax,color=colors,legend=False,width=.8)
                ax.set(ylabel='Estimated contribution (%)',title=label,ylim=(0,100),xlabel='')
            handles,labels=axes[0].get_legend_handles_labels()
            fig.legend(handles,labels,loc='center left',bbox_to_anchor=(.82,.5),fontsize=6)
            warning = 'Exploratory: megakaryocyte reference mismatch unresolved; these estimates do not establish method accuracy.'
            fig.text(.01,.005,warning,fontsize=8)
            fig.tight_layout(rect=(0,.03,.82,1))
            fig.savefig(root/'compositions.pdf',bbox_inches='tight')
            plt.close(fig)
            fig,ax=plt.subplots(figsize=(17,7))
            lim=max(float(diff.abs().to_numpy().max()),1e-9)
            im=ax.imshow(diff.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim)
            ax.set_xticks(range(len(cells)),cells,rotation=90,fontsize=8)
            ax.set_yticks(range(len(samples)),samples,fontsize=8)
            ax.set(title='Candidate CELFIE minus UXM: descriptive differences')
            fig.colorbar(im,ax=ax,label='Percentage points')
            fig.text(.01,.005,warning,fontsize=8)
            fig.tight_layout(rect=(0,.03,1,1))
            fig.savefig(root/'method_differences.pdf',bbox_inches='tight')
            plt.close(fig)
        write_json(root/'interpretation.json', {
            'completed_utc':utc(), 'atlas_sha256':ATLAS_SHA, 'samples':samples,'cell_types':cells,
            'difference_definition':'100 * (CELFIE fraction - preserved UXM fraction)',
            'discrepancy_threshold':None,
            'scope':'Descriptive agreement; no known-truth patient accuracy assessment',
            'limits':[
                'The shared regions were selected for UXM; this is not a neutral comparison of independently optimized marker panels.',
                'UXM fragment patterns and CELFIE CpG counts measure different summaries of the same methylation data.',
                'The reconstructed megakaryocyte counts reproduce checked public donor fingerprints but disagree with the shipped UXM fragment profile at checked regions; this is not an exact shared-reference comparison.',
                'The low-depth megakaryocyte reference can move substantially during native joint fitting; convergence does not establish a credible cell composition.',
                'Reference donors, coverage and native pseudocounts can influence differences.',
                'PB and BM are diagnostic material; these estimates are not validated blast or tumour fractions.',
                'Patient joint fitting can change estimates relative to the initial pilot fit.']})
        finish(root, methods=['uxm','celfie'],samples=samples,cell_types=len(cells))
    except Exception as e:
        finish(root,'failed',error=str(e))
        raise
    return root


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--celfie',type=Path,required=True,help='Completed patient CELFIE run')
    p.add_argument('--fractions',type=Path,required=True,help='Validated sample-by-cell TSV')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--no-plots',action='store_true')
    a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        p.error('Run the comparison in a Slurm allocation')
    print(compare(a))


if __name__=='__main__':
    main()
