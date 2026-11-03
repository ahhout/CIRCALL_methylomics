#!/usr/bin/env python3
"""Observe native CelFiE restarts without editing its source or changing its seeds."""
import argparse
import datetime as dt
import json
from pathlib import Path
import runpy
import sys


def observe(script, output, native_args):
    started = dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')
    script = Path(script).resolve(strict=True)
    output = Path(output)
    with output.open('x') as handle:
        json.dump({'started_utc': started, 'status': 'running'}, handle)
    restarts = []

    def profile(frame, event, value):
        if event != 'return' or frame.f_code.co_name != 'em':
            return
        if Path(frame.f_code.co_filename).resolve() != script:
            return
        if not isinstance(value, tuple) or len(value) != 3:
            return  # Preserve an original native exception without masking it.
        local = frame.f_locals
        alpha, _gamma, likelihood = value
        difference = float(local['alpha_diff'] + local['gamma_diff'])
        restarts.append(dict(seed=int(local['seed']), iterations=int(local['i']) + 1,
                             converged=difference < float(local['convergence_criteria']),
                             final_relative_difference=difference,
                             log_likelihood=float(likelihood), alpha=alpha.tolist()))

    original_argv = sys.argv
    original_profile = sys.getprofile()
    try:
        sys.argv = [str(script), *native_args]
        sys.setprofile(profile)
        namespace = runpy.run_path(str(script), run_name='__main__')
    finally:
        sys.setprofile(original_profile)
        sys.argv = original_argv
    if len(restarts) != namespace['args'].random_restarts:
        raise ValueError('Did not observe every native random restart')
    samples, cells = namespace['samples'], namespace['tissues']
    best = max(range(len(restarts)), key=lambda i: restarts[i]['log_likelihood'])
    # Each range describes variation across the already configured native starts.
    ranges = {}
    for s, sample in enumerate(samples):
        ranges[sample] = {cell: dict(minimum=min(r['alpha'][s][c] for r in restarts),
                                   maximum=max(r['alpha'][s][c] for r in restarts))
                          for c, cell in enumerate(cells)}
    document = dict(scope='Observation of unchanged pinned native CelFiE execution',
                    started_utc=started,
                    completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
                    initial_seed=42, samples=samples, cells=cells, best_restart=best,
                    restarts=restarts, fraction_ranges=ranges,
                    all_starts_converged=all(r['converged'] for r in restarts),
                    best_start_converged=restarts[best]['converged'])
    Path(output).write_text(json.dumps(document, indent=2, allow_nan=False) + '\n')
    return document


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('native_script', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('native_arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from common import output_guard, PROJECT
    output_guard.require_clean_source(PROJECT)
    output_guard.canonical_parent(args.output.absolute().parent)
    output_guard.canonical_parent(Path.cwd())
    observe(args.native_script, args.output, args.native_arguments)


if __name__ == '__main__':
    main()
