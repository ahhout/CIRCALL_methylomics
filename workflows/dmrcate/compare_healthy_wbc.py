#!/usr/bin/env python3
"""Describe the two retained DMRs using exact CpGs and external healthy WBCs."""
import argparse
import bisect
import concurrent.futures
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from output_contract import output_guard, identity


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def table(path, rows, fields=None):
    with Path(path).open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def scan_pat(sample, targets):
    """Stream to EOF for gzip CRC; count all patterns overlapping exact CpGs."""
    started = utc()
    path = Path(sample['pat'])
    if path.stat().st_size != sample['pat_bytes'] or sha(path) != sample['pat_sha256']:
        raise ValueError('PAT differs from the validated download: ' + str(path))
    by_chr = {}
    for i, row in enumerate(targets):
        by_chr.setdefault(row['chr'].encode(), []).append((int(row['cpg_index_1based']), i))
    for values in by_chr.values():
        values.sort()
    indices = {key: [x[0] for x in values] for key, values in by_chr.items()}
    prefixes = tuple(key + b'\t' for key in by_chr)
    methylated = [0] * len(targets)
    coverage = [0] * len(targets)
    records = contributing = 0
    with gzip.open(path, 'rb') as handle:
        for line in handle:
            records += 1
            if not line.startswith(prefixes):
                continue
            fields = line.rstrip(b'\r\n').split(b'\t')
            if len(fields) < 4:
                raise ValueError('Malformed target-chromosome PAT record')
            chrom, start, pattern, count = fields[:4]
            start, count = int(start), int(count)
            if start < 1 or count < 1 or not pattern:
                raise ValueError('Invalid PAT index, multiplicity or pattern')
            positions = indices[chrom]
            left = bisect.bisect_left(positions, start)
            right = bisect.bisect_left(positions, start + len(pattern))
            if left == right:
                continue
            if not re.fullmatch(b'[CT.]+', pattern):
                raise ValueError('Invalid methylation pattern')
            observed = False
            for j in range(left, right):
                position, i = by_chr[chrom][j]
                state = pattern[position - start]
                if state in (67, 84):  # C or T; '.' remains unobserved.
                    coverage[i] += count
                    methylated[i] += count if state == 67 else 0
                    observed = True
            contributing += int(observed)
    return dict(sample=sample['sample'],methylated=methylated,coverage=coverage,
                started_utc=started,completed_utc=utc(),gzip_CRC='passed: streamed to EOF',
                sha256=sample['pat_sha256'],records_streamed=records,
                records_contributing_to_targets=contributing)


def summarize(rows):
    groups = {}
    for row in rows:
        groups.setdefault((row['region'],row['sample']), []).append(row)
    result = []
    for (region, sample), values in groups.items():
        cov = [int(x['coverage']) for x in values]
        complete = all(x > 0 for x in cov)
        mean = statistics.mean(int(x['methylated'])/int(x['coverage']) for x in values) if complete else None
        result.append(dict(region=region,sample=sample,group=values[0]['group'],
            patient_or_donor=values[0]['patient_or_donor'],assay=values[0]['assay'],
            target_cpgs=len(values),observed_cpgs=sum(x>0 for x in cov),
            mean_fraction=mean if mean is not None else 'NA',
            minimum_coverage=min(cov),total_coverage=sum(cov),
            all_cpgs_at_least_3x=all(x>=3 for x in cov)))
    return result


def run(args):
    started = utc()
    out = args.output.absolute()
    if not re.fullmatch(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-[A-Za-z0-9][A-Za-z0-9._-]*', out.name):
        raise ValueError('Output requires a UTC date-first run name')
    dt.datetime.strptime(out.name[:16], '%Y-%m-%d_%H-%M')
    if not 1 <= args.threads <= 3:
        raise ValueError('Use one to three PAT workers, within the eight-CPU limit')
    if out.exists() or out.is_symlink():
        raise ValueError('Refusing to overwrite an existing comparison')
    patient = args.patient_run.resolve()
    controls = args.controls.resolve()
    staged = json.loads((controls/'STAGED.json').read_text())
    cfg = json.loads((patient/'config.json').read_text())
    if args.runtime.resolve() != Path(cfg['runtime']).resolve():
        raise ValueError('Use the patient run\'s frozen runtime')
    for executable in (args.tabix,args.apptainer):
        if not executable.is_file() or not os.access(executable,os.X_OK):
            raise ValueError('Required executable is unavailable: '+str(executable))
    subprocess.run(['sha256sum','-c',str(args.runtime/'sif.sha256')],check=True)
    if staged['status'] != 'passed' or staged['build'] != 'hg38' or not (patient/'COMPLETE').is_file():
        raise ValueError('Require validated hg38 controls and a completed patient run')
    if len(staged['samples']) != 3:
        raise ValueError('This comparison expects the three selected healthy WBC donors')
    dictionary = Path(cfg['namespace']).resolve()
    fingerprint = json.loads((patient/'input_fingerprint.json').read_text())
    dictionary_sha = sha(dictionary)
    if dictionary_sha != staged['dictionary_sha256'] or fingerprint[str(dictionary)] != dictionary_sha:
        raise ValueError('Control and patient CpG dictionaries differ')
    targets = []
    regions = []
    for key in ('cfDNA_minus_BM','cfDNA_minus_PB'):
        selected = read_table(patient/'comparisons'/key/'regions.tsv')
        if len(selected) != 1:
            raise ValueError('Expected one retained DMR per cfDNA contrast')
        row = selected[0]
        expected = {'cfDNA_minus_BM':'chr1:629824-629882',
                    'cfDNA_minus_PB':'chr20:61252799-61252885'}
        if row['coord'] != expected[key]:
            raise ValueError('Candidate regions differ from the explicitly selected follow-up')
        chrom, interval = row['coord'].split(':')
        start, end = interval.split('-')
        text = subprocess.check_output([str(args.tabix),str(dictionary),f'{chrom}:{start}-{end}'],text=True)
        sites = [line.split() for line in text.splitlines()]
        if len(sites) != int(row['no.cpgs']):
            raise ValueError('Dictionary interval differs from the retained tested CpG set')
        for ch, pos, index in sites:
            targets.append(dict(region=row['coord'],source_contrast=key,chr=ch,
                                position_1based=int(pos),cpg_index_1based=int(index)))
        regions.append(row['coord'])
    output_guard.reserve_path(out, identity(Path(__file__), [patient/'config.json', controls/'STAGED.json']))
    (out/'code').mkdir(); (out/'tmp').mkdir()
    for name in ('compare_healthy_wbc.py','extract_candidate_counts.R'):
        shutil.copy2(Path(__file__).with_name(name),out/'code'/name)
    table(out/'targets.tsv', targets)
    manifest = dict(started_utc=started,status='running',patient_run=str(patient),
        controls=str(controls),argv=sys.argv,node=os.uname().nodename,
        execution='Direct current-node check; three single-thread PAT workers and one serial R extraction.',
        dictionary_sha256=dictionary_sha,genome='hg38',target_regions=regions,
        decision='Exact existing target CpGs; no new tests, pooled fitting, filtering or imputation.',
        mean_rule='Arithmetic mean of per-CpG fractions; NA if any predefined CpG is uncovered. No extra coverage filter.',
        healthy_plasma_specificity='On hold; these are external healthy adult WBC references.',
        input_sha256={str(p):sha(p) for p in [patient/'config.json',patient/'prepared/samples.tsv',
            controls/'STAGED.json',controls/'code/healthy_blood_proposal.json',
            args.runtime/'sif.sha256',args.runtime/'packages.lock.tsv']},
        code_sha256={p.name:sha(p) for p in (out/'code').iterdir()})
    def save_manifest():
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    save_manifest()
    try:
        env = dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        command = ['/usr/bin/time','-v','-o',str(out/'patient-resources.txt'),
            str(args.apptainer),'exec','--cleanenv','--no-home','--containall',
            '--bind',f'{patient}:{patient}:ro','--bind',f'{args.runtime}:{args.runtime}:ro',
            '--bind',f'{out}:{out}:rw','--env',f'TMPDIR={out}/tmp,OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,TZ=UTC',
            str(args.runtime/'dmrcate.sif'),'Rscript',str(out/'code/extract_candidate_counts.R'),
            str(patient),str(out/'targets.tsv'),str(out)]
        with (out/'patient-extraction.log').open('x') as log:
            subprocess.run(command,check=True,env=env,stdout=log,stderr=subprocess.STDOUT)
        samples = {x['sample']:x for x in read_table(patient/'prepared/samples.tsv')}
        rows = read_table(out/'patient_cpg_counts.tsv')
        for row in rows:
            row.update(group=samples[row['sample']]['compartment'],
                patient_or_donor=samples[row['sample']]['patient'],assay='EM-seq')
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as pool:
            futures = [pool.submit(scan_pat,sample,targets) for sample in staged['samples']]
            control_results = [future.result() for future in futures]
        for sample, result in zip(staged['samples'],control_results):
            for i,target in enumerate(targets):
                rows.append(dict(target,sample=sample['sample'],methylated=result['methylated'][i],
                    coverage=result['coverage'][i],group='Healthy WBC',
                    patient_or_donor=sample['donor_label'],assay='WGBS'))
        table(out/'cpg_counts.tsv',rows)
        summaries = summarize(rows)
        # Confirm the patient extraction reproduces every original reported mean.
        for key in ('cfDNA_minus_BM','cfDNA_minus_PB'):
            for original in read_table(patient/'comparisons'/key/'regional_sample_methylation.tsv'):
                row = next(x for x in summaries if x['region']==original['region'] and x['sample']==original['sample'])
                if abs(float(row['mean_fraction'])-float(original['mean_fraction'])) > 1e-12:
                    raise ValueError('Patient methylation differs from the original result')
        table(out/'region_summary.tsv',summaries)
        os.environ['MPLCONFIGDIR'] = str(out/'tmp/matplotlib')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1,2,figsize=(10,4.8),sharey=True)
        groups = ['Healthy WBC','PB','BM','cfDNA']
        people = sorted({x['patient_or_donor'] for x in summaries if x['group']!='Healthy WBC'})
        for ax,region,label in zip(axes,regions,['MTND2P28 region (3 CpGs)','CDH4 region (12 CpGs)']):
            values = [x for x in summaries if x['region']==region]
            for i,person in enumerate(people):
                subset = [next(x for x in values if x['patient_or_donor']==person and x['group']==g) for g in groups[1:]]
                y = [100*float(x['mean_fraction']) if x['mean_fraction']!='NA' else float('nan') for x in subset]
                ax.plot([1,2,3],y,'o-',label=person.split('-')[-1],color=f'C{i}',linewidth=1.2)
            healthy = [x for x in values if x['group']=='Healthy WBC']
            for dx,row in zip([-.08,0,.08],healthy):
                if row['mean_fraction']!='NA': ax.scatter(dx,100*float(row['mean_fraction']),color='black',marker='s',zorder=4)
            ax.set_xticks(range(4),['Healthy\nWBC','Patient\nPB','Patient\nBM','Patient\ncfDNA'])
            ax.set_title(label+'\n'+region,fontsize=10); ax.set_ylim(-3,103)
            ax.grid(axis='y',alpha=.25)
        axes[0].set_ylabel('Mean CpG methylation (%)'); axes[1].legend(title='Patient',fontsize=8)
        fig.suptitle('External WGBS controls and patient EM-seq: descriptive comparison',fontsize=11)
        fig.tight_layout(); fig.savefig(out/'candidate_methylation.png',dpi=180);fig.savefig(out/'candidate_methylation.pdf')
        plt.close(fig)
        manifest.update(status='completed',completed_utc=utc(),control_integrity=control_results,
            patient_means_match_original=True,
            missing_region_means=sum(x['mean_fraction']=='NA' for x in summaries),
            limitation='Adult female WBC WGBS controls are not matched healthy plasma; chemistry, processing and age/sex differ.')
        save_manifest(); (out/'COMPLETE').write_text(manifest['completed_utc']+'\n')
        output_guard.finish(out, True)
        print(json.dumps(dict(output=str(out),status='completed',regions=len(regions),
                              controls=len(control_results),region_summaries=len(summaries))))
    except Exception as error:
        manifest.update(status='failed',failed_utc=utc(),error=str(error));save_manifest()
        output_guard.finish(out, False)
        raise


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--patient-run',type=Path,required=True)
    parser.add_argument('--controls',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--runtime',type=Path,required=True)
    parser.add_argument('--tabix',type=Path,required=True)
    parser.add_argument('--apptainer',type=Path,required=True)
    parser.add_argument('--threads',type=int,default=3)
    run(parser.parse_args())
