#!/usr/bin/env python3
"""Generate two known mixtures and exercise native routes, missingness and resume."""
import argparse
import sys
import gzip
import os
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import *
from native import verify_runtime, invoke, native_fractions


def fixture_inputs(root):
    markers, _ = atlas_rows()
    markers = [m for m in markers if m['chr'] == 'chr1' and int(m['endCpG']) - int(m['startCpG']) >= 5][:20]
    if len(markers) != 20:
        raise ValueError('Insufficient fixture intervals')
    fixture = root / 'inputs'
    fixture.mkdir()
    cells = ['fixture_A', 'fixture_B']
    u_atlas, beta_atlas = [], []
    for i, marker in enumerate(markers):
        # Non-complementary profiles avoid an intercept-collinear RPC design.
        a, b = [(.9, .1), (.2, .8), (.8, .3), (.4, .9)][i % 4]
        u_atlas.append(dict(list(marker.items())[:8], target=cells[i % 2], fixture_A=a, fixture_B=b))
        # Only this toy model has exclusively fully methylated/unmethylated reads,
        # so beta = 1-U here. The real reference builder aggregates PAT counts.
        beta_atlas.append(dict(chr=marker['chr'], start=marker['start'], end=marker['end'], fixture_A=1-a, fixture_B=1-b))
    # Native UXM infers its minimum read length from the atlas's .l4. token.
    table(fixture / 'Atlas.U25.l4.hg38.fixture.tsv', u_atlas, list(markers[0])[:8] + cells)
    table(fixture / 'atlas.beta.tsv', beta_atlas, ['chr', 'start', 'end'] + cells)
    sequences = []
    for m in markers:
        region = f'{m["chr"]}:{int(m["start"])+1}-{m["end"]}'
        sequence = run(['samtools', 'faidx', GENOME / 'genome.fa.gz', region], capture_output=True, text=True).stdout
        sequences.append(''.join(sequence.splitlines()[1:]).upper())
    bam_sheet, cov_sheet = [], []
    chroms = [line.split() for line in (GENOME / 'chrome.size').read_text().splitlines()]
    for sample, mixture in [('fixture_a', 7), ('fixture_b', 3)]:
        sam = fixture / (sample + '.sam')
        cov = fixture / (sample + '.cov.gz')
        with sam.open('w') as handle, gzip.open(cov, 'wt') as coverage:
            handle.write('@HD\tVN:1.6\tSO:unsorted\n')
            for chrom, size in chroms:
                handle.write(f'@SQ\tSN:{chrom}\tLN:{size}\n')
            for i, (m, sequence) in enumerate(zip(markers, sequences)):
                a, b = float(u_atlas[i]['fixture_A']), float(u_atlas[i]['fixture_B'])
                u_reads = round(100 * (a * mixture / 10 + b * (1-mixture / 10)))
                cpg = [j for j in range(len(sequence)-1) if sequence[j:j+2] == 'CG']
                if len(cpg) < 4:
                    raise ValueError('Fixture read does not cover four CpGs')
                for k in range(100):
                    methylated = k >= u_reads
                    xm = ['.'] * len(sequence)
                    bases = list(sequence)
                    for j in cpg:
                        xm[j] = 'Z' if methylated else 'z'
                        if not methylated:
                            bases[j] = 'T'
                    handle.write(f'{sample}_{i}_{k}\t0\t{m["chr"]}\t{int(m["start"])+1}\t60\t{len(sequence)}M\t*\t0\t0\t{"".join(bases)}\t{"I"*len(sequence)}\tXM:Z:{"".join(xm)}\tXR:Z:CT\tXG:Z:CT\n')
                # One missing region tests native per-sample missingness handling.
                if sample == 'fixture_b' and i == len(markers)-1:
                    continue
                for j in cpg:
                    position = int(m['start']) + 1 + j
                    coverage.write(f'{m["chr"]}\t{position}\t{position}\t{100-u_reads}\t{100-u_reads}\t{u_reads}\n')
        bam = fixture / (sample + '.bam')
        run(['samtools', 'sort', '-o', bam, sam], root / 'commands.jsonl')
        run(['samtools', 'index', bam], root / 'commands.jsonl')
        bam_sheet.append(dict(name=sample, bam=str(bam), bai=str(bam)+'.bai'))
        cov_sheet.append(dict(name=sample, sample=str(cov)))
    table(fixture / 'bams.csv', bam_sheet, ['name', 'bam', 'bai'], ',')
    table(fixture / 'coverage.csv', cov_sheet, ['name', 'sample'], ',')
    return ['fixture_a', 'fixture_b'], cells


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runtime', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Use validate_runtime.sh through sbatch')
    root = new_run(args.output)
    try:
        runtime = verify_runtime(args.runtime)
        samples, cells = fixture_inputs(root)
        checks = []
        # Exercise both methods independently even if one has an upstream defect.
        errors = []
        for method, sheet, atlas in [('uxm', 'bams.csv', 'Atlas.U25.l4.hg38.fixture.tsv'), ('epidish', 'coverage.csv', 'atlas.beta.tsv')]:
            try:
                result = invoke(runtime, root, method, root / 'inputs' / sheet, root / 'inputs' / atlas, 'fixture', fixtures=True)
                fractions = native_fractions(result, method, samples, cells)
                if not np.allclose(fractions['fixture_A'], [.7, .3], atol=.02, rtol=0):
                    raise ValueError(f'{method}: synthetic mixture estimates differ from expected fractions')
                before = fractions.copy()
                # A distinct report destination exercises cache reuse, preserving first-attempt traces.
                resumed = invoke(runtime, root, method, root / 'inputs' / sheet, root / 'inputs' / atlas, 'recovery', resume=True, fixtures=True)
                after = native_fractions(resumed, method, samples, cells)
                if not before.equals(after):
                    raise ValueError('Resume changed native estimates')
                trace = rows(resumed / 'pipeline_info/execution_trace_recovery.txt')
                if any(r['status'] != 'CACHED' for r in trace):
                    raise ValueError('Resume did not reuse all completed native processes')
                checks.append({'method': method, 'status': 'passed', 'checks': ['known mixtures', 'sample labels', 'fraction bounds', 'resume']})
            except Exception as exc:
                errors.append({'method': method, 'error': str(exc)})
        write_json(root / 'checks.json', {'passed': checks, 'failed': errors})
        if errors:
            raise ValueError('Native fixture failure: ' + json.dumps(errors))
        write_json(root / 'VALIDATED.json', {'runtime_lock': sha(runtime / 'runtime.lock.json'), 'completed_utc': utc(),
                   'checks': checks, 'execution': 'Slurm; NXF_OFFLINE=true; pre-staged plugin/images; invalid HTTP proxies'})
        finish(root)
    except Exception as exc:
        finish(root, 'failed', error=str(exc))
        raise


if __name__ == '__main__':
    main()
