#!/usr/bin/env python3
"""Validate native CelFiE with synthetic counts only; no patient reference is inferred."""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from common import *
from native import checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--runtime', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Submit run_celfie_fixture.sh through sbatch')
    root = new_run(args.output)
    try:
        runtime = args.runtime.resolve(strict=True)
        source = runtime / ('DecoNFlow-' + COMMIT)
        lock = read_json(runtime / 'runtime.lock.json')
        required = ['images/python.sif'] + [
            str((source / relative).relative_to(runtime))
            for relative in ('bin/celfie.py', 'modules/celfie/deconv/main.nf')]
        for relative in required:
            if sha(runtime / relative) != lock[relative]:
                raise ValueError(f'Pinned CelFiE dependency changed: {relative}')
        fixture = root / 'counts.tsv'
        columns = ['CHR', 'START', 'END', 'fixture_a_meth', 'fixture_a_depth',
                   'fixture_b_meth', 'fixture_b_depth', 'REF_CHR', 'REF_START',
                   'REF_END', 'A_meth', 'A_depth', 'B_meth', 'B_depth']
        records = []
        for i in range(40):
            a, b = (.1, .9) if i % 2 else (.8, .2)
            records.append(['chr1', i * 100, i * 100 + 50,
                            round(10000 * (.7 * a + .3 * b)), 10000,
                            round(10000 * (.3 * a + .7 * b)), 10000,
                            'chr1', i * 100, i * 100 + 50,
                            round(a * 100000), 100000, round(b * 100000), 100000])
        pd.DataFrame(records, columns=columns).to_csv(fixture, sep='\t', index=False)
        parameters = dict(deconflow_source=str(source), count_table=str(fixture),
                          outdir=str(root / 'results'), celfie_maxiter=1000, unknown=0,
                          converg=0.0001, celfie_randrest=10, parall_job=1)
        write_json(root / 'params.json', parameters)
        config = root / 'fixture.config'
        config.write_text(f"""
process {{
    executor = 'slurm'
    queue = 'shortq'
    cpus = 2
    memory = '4 GB'
    time = '20m'
    errorStrategy = 'terminate'
    withName: RUN_DECONV {{ container = '{runtime}/images/python.sif' }}
    beforeScript = 'mkdir -p .tmp'
    containerOptions = {{ '--bind ' + task.workDir + '/.tmp:/tmp --bind {runtime}:{runtime}:ro --env TMPDIR=/tmp' }}
    publishDir = [path: '{root}/results', mode: 'copy']
}}
apptainer {{ enabled = true; autoMounts = true; cacheDir = '/path/to/hpc-scratch/nextflow/apptainer_cache' }}
env {{
    PATH = '{source}/bin:/usr/local/bin:/usr/bin:/bin'
    TMPDIR = '/tmp'
    OMP_NUM_THREADS = '1'
    OPENBLAS_NUM_THREADS = '1'
    http_proxy = 'http://127.0.0.1:9'
    https_proxy = 'http://127.0.0.1:9'
}}
trace {{
    enabled = true
    file = params.trace_path
    fields = 'task_id,hash,native_id,name,status,exit,duration,peak_rss,workdir'
}}
cleanup = false
""")
        env = dict(os.environ, NXF_OFFLINE='true', NXF_ANSI_LOG='false',
                   NXF_PLUGINS_DIR=str(runtime / 'plugins'))
        results = {}
        for phase in ('first', 'recovery'):
            command = ['nextflow', 'run', WORKFLOW / 'tools/celfie_fixture.nf',
                       '-c', config, '-params-file', root / 'params.json',
                       '-work-dir', WORK, '--trace_path', root / f'{phase}.trace.tsv']
            if phase == 'recovery':
                command.extend(['-resume', resume_nextflow_session(root)])
            with (root / f'{phase}.log').open('w') as log:
                run(command, root / 'commands.jsonl', cwd=root, env=env,
                    stdout=log, stderr=subprocess.STDOUT)
            save_nextflow_session(root)
            output = root / 'results/1_tissue_proportions.txt'
            fractions = checked(pd.read_csv(output, sep='\t', index_col=0),
                                ['fixture_a', 'fixture_b'], ['A', 'B'])
            if not np.allclose(fractions['A'], [.7, .3], atol=.02, rtol=0):
                raise ValueError('Native CelFiE did not recover the known mixtures')
            trace = rows(root / f'{phase}.trace.tsv')
            wanted = 'CACHED' if phase == 'recovery' else 'COMPLETED'
            if len(trace) != 1 or trace[0]['status'] != wanted:
                raise ValueError(f'Unexpected {phase} trace: {trace}')
            results[phase] = sha(output)
        if results['first'] != results['recovery']:
            raise ValueError('Cached recovery changed fractions')
        write_json(root / 'VALIDATED.json',
                   dict(completed_utc=utc(), scope='synthetic native CelFiE and cached recovery only',
                        dependencies={p: lock[p] for p in required}, outputs=results,
                        patient_reference_validated=False))
        finish(root)
    except Exception as exc:
        finish(root, 'failed', error=str(exc))
        raise


if __name__ == '__main__':
    main()
