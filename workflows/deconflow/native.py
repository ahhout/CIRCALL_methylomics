"""Native Nextflow invocation and runtime verification shared by fixtures and analysis."""
import os
import numpy as np
import pandas as pd
from common import *


def verify_runtime(runtime):
    runtime = Path(runtime).resolve(strict=True)
    if not (runtime / 'STAGED').is_file():
        raise ValueError('Runtime staging did not complete')
    for relative, expected in read_json(runtime / 'runtime.lock.json').items():
        if sha(runtime / relative) != expected:
            raise ValueError(f'Frozen runtime changed: {relative}')
    for alias in ('hg38', 'default'):
        if (runtime / 'references' / alias).resolve(strict=True) != GENOME:
            raise ValueError(f'Isolated reference alias changed: {alias}')
    if sha(GENOME / 'CpG.bed.gz') != CPG_SHA:
        raise ValueError('Canonical hg38 CpG dictionary changed')
    return runtime


def configuration(runtime, target, fixtures=False):
    # Snapshot explicit configuration instead of generating Groovy in Python.
    if any(any(c in str(p) for c in "'\n,") for p in (runtime, PROJECT, GENOME, WORK)):
        raise ValueError('Unsupported characters in runtime paths')
    config = (WORKFLOW / 'config/flamingo.config').read_bytes()
    if fixtures:
        config += b'\n' + (WORKFLOW / 'config/fixtures.config').read_bytes()
    target.write_bytes(config)


def invoke(runtime, run_root, method, sample_sheet, atlas, phase, resume=False, fixtures=False,
           fragments=False, result_override=None):
    if fragments and method not in ('uxm', 'counts'):
        raise ValueError('The local PAT entry point supports UXM and its counter only; EpiDISH is retired')
    require_incomplete(run_root)
    result = result_override or run_root / phase / method
    if not Path(result).resolve().is_relative_to(Path(run_root).resolve()):
        raise ValueError('Native publication must stay inside its reserved run bundle')
    launch = run_root / 'native' / method
    if fragments:
        launch = launch / phase
    launch.mkdir(parents=True, exist_ok=True)
    config = run_root / 'flamingo.config'
    if not config.exists():
        configuration(runtime, config, fixtures)
    result = result_override or run_root / phase / method
    params = read_json(WORKFLOW / 'config' / f'params.{method}.json')
    params.update(outdir=str(result), trace_report_suffix=phase)
    params['test_bams' if method == 'uxm' else 'test_set'] = str(sample_sheet)
    params['uxm_atlas' if method == 'uxm' else 'ref_matrix'] = str(atlas)
    if fragments:
        params.update(fragment_method=method, fragment_samples=str(sample_sheet),
                      fragment_reference=str(atlas), deconflow_source=str(runtime / ('DecoNFlow-' + COMMIT)))
    params_path = launch / (phase + '.params.json')
    if params_path.exists() and read_json(params_path) != params:
        raise ValueError('Refusing changed native parameters on resume')
    write_json(params_path, params)
    env = dict(os.environ, NXF_OFFLINE='true', NXF_PLUGINS_DIR=str(runtime / 'plugins'),
               NXF_PLUGINS_ALLOWED='nf-schema', NXF_ANSI_LOG='false', NXF_WORK=str(WORK),
               NXF_ENABLE_VIRTUAL_THREADS='false',
               DECONFLOW_RUNTIME=str(runtime), DECONFLOW_GENOME=str(GENOME))
    entrypoint = WORKFLOW / 'fragments.nf' if fragments else runtime / ('DecoNFlow-' + COMMIT) / 'main.nf'
    command = ['nextflow', 'run', entrypoint,
               '-c', config, '-params-file', params_path,
               '-work-dir', WORK, '-name', f'deconflow_{method}_{phase}_{os.environ.get("SLURM_JOB_ID", "test")}']
    if fragments:
        command.extend(['-c', WORKFLOW / 'config/fragments.config'])
    else:
        command.extend(['-profile', 'apptainer'])
    if resume:
        command.extend(['-resume', resume_nextflow_session(launch)])
    log = launch / f'{phase}-{os.environ.get("SLURM_JOB_ID", "test")}.log'
    with log.open('a') as handle:
        try:
            run(command, run_root / 'commands.jsonl', cwd=launch, env=env, stdout=handle, stderr=subprocess.STDOUT)
        except BaseException:
            if (launch / '.nextflow.log').is_file():
                save_nextflow_session(launch)
            raise
    save_nextflow_session(launch)
    trace = result / 'pipeline_info' / f'execution_trace_{phase}.txt'
    if not trace.exists():
        raise ValueError('Native workflow produced no execution trace')
    entries = rows(trace)
    if not entries or any(r.get('status') not in ('COMPLETED', 'CACHED') for r in entries):
        raise ValueError('Native trace contains incomplete tasks')
    return result


def checked(frame, samples, cells):
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError('Duplicate sample/cell labels')
    if set(frame.index) != set(samples) or set(frame.columns) != set(cells):
        raise ValueError('Missing/extra sample or cell labels in fraction table')
    frame = frame.loc[samples, cells].apply(pd.to_numeric, errors='raise')
    if not np.isfinite(frame.values).all() or (frame.values < -1e-8).any() or (frame.values > 1 + 1e-8).any():
        raise ValueError('Nonfinite or out-of-range fractions')
    # Numerical validation only, not a scientific discrepancy threshold.
    if not np.allclose(frame.sum(axis=1), 1, atol=1e-5, rtol=0):
        raise ValueError('Fractions do not sum to one within output rounding precision')
    return frame


def native_fractions(root, method, samples, cells):
    relative = 'uxm/output.csv' if method == 'uxm' else 'epidish/epidish_res_RPC.csv'
    path = root / 'deconvolution_results/reference_based' / relative
    frame = pd.read_csv(path, index_col=0)
    return checked(frame, samples, cells)
