#!/usr/bin/env python3
"""Observe unchanged native CelFiE and compare terminal restart scores diagnostically.

This helper never substitutes a restart into the native outputs. Observed-data
scores use the reference counts after native pseudocount handling, with their
contribution included once. No biological or convergence acceptance is granted.
"""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import runpy
import sys

import numpy as np


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def identity(path):
    path = Path(path).resolve(strict=True)
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def _counts(methylated, depths, label):
    methylated = np.asarray(methylated, dtype=float)
    depths = np.asarray(depths, dtype=float)
    if (methylated.ndim != 2 or methylated.shape != depths.shape
            or not np.isfinite(methylated).all() or not np.isfinite(depths).all()
            or (methylated < 0).any() or (depths < methylated).any()):
        raise ValueError(f'Invalid {label} methylated/depth arrays')
    return methylated, depths


def _weighted_log(counts, probability):
    """Use 0 log(0) = 0; positive count at probability zero is impossible."""
    positive = counts > 0
    with np.errstate(divide='ignore'):
        return float(np.sum(counts[positive] * np.log(probability[positive])))


def observed_log_likelihood(x, x_depths, y, y_depths, alpha, gamma):
    """Return total, sample and reference scores, omitting binomial constants."""
    x, x_depths = _counts(x, x_depths, 'sample')
    y, y_depths = _counts(y, y_depths, 'reference')
    alpha = np.asarray(alpha, dtype=float)
    gamma = np.asarray(gamma, dtype=float)
    if (alpha.shape != (x.shape[0], y.shape[0]) or gamma.shape != y.shape
            or x.shape[1] != y.shape[1] or not np.isfinite(alpha).all()
            or not np.isfinite(gamma).all() or (alpha < 0).any()
            or (alpha > 1).any() or (gamma < 0).any() or (gamma > 1).any()
            or not np.allclose(alpha.sum(axis=1), 1, rtol=0, atol=1e-10)):
        raise ValueError('Invalid terminal alpha/gamma arrays')
    probability = alpha @ gamma
    if (probability < 0).any() or (probability > 1).any():
        raise ValueError('Mixture methylation probability outside [0, 1]')
    sample = (_weighted_log(x, probability)
              + _weighted_log(x_depths - x, 1 - probability))
    reference = (_weighted_log(y, gamma)
                 + _weighted_log(y_depths - y, 1 - gamma))
    return dict(total=sample + reference, sample=sample, reference=reference)


def _json_score(value):
    if np.isnan(value) or np.isposinf(value):
        raise ValueError('Invalid restart score')
    return None if np.isneginf(value) else float(value)


def _run_root(output):
    for parent in output.parents:
        if re.fullmatch(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-.+-v[1-9]\d*', parent.name):
            dt.datetime.strptime(parent.name[:16], '%Y-%m-%d_%H-%M')
            return parent
    raise ValueError('Diagnostic output needs an existing YYYY-MM-DD_HH-MM-<purpose>-vN run ancestor')


def observe(script, output, native_args):
    """Profile native ``em`` returns without changing arrays, RNG or outputs."""
    started = utc()
    script = Path(script).resolve(strict=True)
    output = Path(output).absolute()
    if not output.parent.is_dir():
        raise ValueError('Create the dated run and diagnostic directory before invoking the observer')
    run_root = _run_root(output)
    states_path = output.with_name(output.stem + '.states.npz')
    if output.exists() or output.is_symlink() or states_path.exists() or states_path.is_symlink():
        raise ValueError('Refusing to overwrite existing score diagnostics or terminal states')
    native_args = [str(item) for item in native_args]
    if len(native_args) < 3:
        raise ValueError('Native arguments must include input path, output directory and sample count')
    native_options = argparse.ArgumentParser(add_help=False)
    native_options.add_argument('-p', '--parallel_job_id', type=int, default=1)
    parsed_options, _ = native_options.parse_known_args(native_args[3:])
    for suffix in ('tissue_proportions', 'methylation_proportions'):
        native_output = Path.cwd() / f'{parsed_options.parallel_job_id}_{suffix}.txt'
        if native_output.exists() or native_output.is_symlink():
            raise ValueError(f'Refusing to overwrite native output; use a fresh task directory: {native_output}')
    document = dict(
        schema_version=1, status='running', started_utc=started,
        scope='Diagnostic observation of unchanged native CelFiE terminal states; no acceptance',
        run_root=str(run_root), native_script=identity(script),
        observer_script=identity(__file__), input=identity(native_args[0]),
        native_argv=[str(script), *native_args], native_outputs_unchanged_by_observer=True,
        score_definition='Observed-data binomial log likelihood without combinatorial constants; reference counted once',
        reference_score_basis='Effective native reference counts after pseudocount handling',
        negative_infinity_encoding='A null score with negative_infinity=true means an impossible observation',
        biological_accuracy_validated=False, acceptance_granted=False)
    # Exclusive creation records execution start before any native work begins.
    with output.open('x') as handle:
        json.dump(document, handle, indent=2, allow_nan=False)
        handle.write('\n')
    captured = []
    callback_errors = []

    def profile(frame, event, value):
        if (event != 'return' or frame.f_code.co_name != 'em'
                or Path(frame.f_code.co_filename).resolve() != script):
            return
        if not isinstance(value, tuple) or len(value) != 3:
            return  # Do not mask a native exception.
        try:
            local = frame.f_locals
            alpha, gamma, native_score = value
            entry = {name: np.array(local[name], copy=True)
                     for name in ('x', 'x_depths', 'y', 'y_depths')}
            entry.update(alpha=np.array(alpha, copy=True), gamma=np.array(gamma, copy=True),
                         candidate_alpha=np.array(local['a'], copy=True),
                         candidate_gamma=np.array(local['g'], copy=True),
                         native_score=float(native_score), seed=int(local['seed']),
                         iterations=int(local['i']) + 1,
                         threshold=float(local['convergence_criteria']),
                         alpha_relative_change=float(local['alpha_diff']),
                         gamma_relative_change=float(local['gamma_diff']))
            captured.append(entry)
        except Exception as exc:
            # Complete the native calculation even if this observer cannot inspect it.
            callback_errors.append(f'{type(exc).__name__}: {exc}')

    original_argv, original_profile = sys.argv, sys.getprofile()
    original_np_errors = np.geterr()
    try:
        try:
            sys.argv = [str(script), *native_args]
            sys.setprofile(profile)
            namespace = runpy.run_path(str(script), run_name='__main__')
        finally:
            sys.setprofile(original_profile)
            sys.argv = original_argv
            np.seterr(**original_np_errors)
        if callback_errors:
            raise ValueError(f'Native restart capture failed: {callback_errors}')
        if not captured or len(captured) != namespace['args'].random_restarts:
            raise ValueError('Did not observe every native random restart')
        if identity(script) != document['native_script'] or identity(native_args[0]) != document['input']:
            raise ValueError('Native source or input changed during execution')
        samples, cells = list(namespace['samples']), list(namespace['tissues'])
        first = captured[0]
        for name, label in (('y', 'reference'), ('y_depths', 'reference'),
                            ('x', 'sample'), ('x_depths', 'sample')):
            if not all(np.array_equal(item[name], first[name]) for item in captured):
                document[f'effective_{label}_counts_identical'] = False
                raise ValueError(f'Effective {label} counts differ across restarts; scores are not comparable')
        document.update(effective_reference_counts_identical=True,
                        effective_sample_counts_identical=True)
        summaries, observed_scores = [], []
        for index, item in enumerate(captured):
            if item['alpha'].shape != (len(samples), len(cells)):
                raise ValueError('Observed terminal fractions do not match native labels')
            scores = observed_log_likelihood(item['x'], item['x_depths'], item['y'],
                                            item['y_depths'], item['alpha'], item['gamma'])
            observed_scores.append(scores['total'])
            difference = item['alpha_relative_change'] + item['gamma_relative_change']
            if (not np.isfinite(difference) or not np.isfinite(item['threshold'])
                    or item['threshold'] <= 0):
                raise ValueError('Nonfinite native convergence diagnostic or invalid threshold')
            summaries.append(dict(
                restart=index, seed=item['seed'], iterations=item['iterations'],
                converged=difference < item['threshold'], convergence_threshold=item['threshold'],
                final_relative_difference=difference,
                alpha_relative_change=item['alpha_relative_change'],
                gamma_relative_change=item['gamma_relative_change'],
                candidate_alpha_max_absolute_difference_from_returned=float(np.max(np.abs(item['candidate_alpha'] - item['alpha']))),
                candidate_gamma_max_absolute_difference_from_returned=float(np.max(np.abs(item['candidate_gamma'] - item['gamma']))),
                native_score=_json_score(item['native_score']),
                native_score_negative_infinity=bool(np.isneginf(item['native_score'])),
                observed_log_likelihood=_json_score(scores['total']),
                observed_log_likelihood_negative_infinity=bool(np.isneginf(scores['total'])),
                sample_log_likelihood=_json_score(scores['sample']),
                reference_log_likelihood=_json_score(scores['reference']),
                alpha=item['alpha'].tolist()))
        native_best = max(range(len(captured)), key=lambda i: captured[i]['native_score'])
        observed_best = max(range(len(captured)), key=lambda i: observed_scores[i])
        converged = [i for i, item in enumerate(summaries) if item['converged']]
        converged_best = max(converged, key=lambda i: observed_scores[i]) if converged else None
        # All arrays are copies of returned/local states; no rescored state is fitted.
        arrays = {name: np.stack([item[name] for item in captured]) for name in
                  ('alpha', 'gamma', 'candidate_alpha', 'candidate_gamma')}
        arrays.update({name: first[name] for name in ('x', 'x_depths', 'y', 'y_depths')})
        arrays.update(seeds=np.array([item['seed'] for item in captured]),
                      native_scores=np.array([item['native_score'] for item in captured]),
                      observed_scores=np.array(observed_scores))
        with states_path.open('xb') as handle:
            np.savez_compressed(handle, **arrays)
        native_outputs = {name: identity(namespace[name]) for name in
                          ('output_alpha_file', 'output_gamma_file') if name in namespace}
        document.update(
            status='complete', completed_utc=utc(), samples=samples, cells=cells,
            initial_seed=namespace.get('initial_seed'), restarts=summaries,
            best_native_restart=native_best, best_observed_restart=observed_best,
            best_converged_observed_restart=converged_best,
            native_best_start_converged=summaries[native_best]['converged'],
            best_observed_start_converged=summaries[observed_best]['converged'],
            all_starts_converged=all(item['converged'] for item in summaries),
            observed_best_is_nonconverged=not summaries[observed_best]['converged'],
            selection_policy='Native output retained; observed and converged rankings are descriptive only; no fallback accepted',
            states=identity(states_path),
            state_array_note='y/y_depths and x/x_depths stored once after equality across all captured restarts was checked',
            native_outputs=native_outputs)
    except BaseException as exc:
        document.update(status='failed', ended_utc=utc(), error=f'{type(exc).__name__}: {exc}',
                        captured_restarts=len(captured))
        raise
    finally:
        output.write_text(json.dumps(document, indent=2, allow_nan=False) + '\n')
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
