"""Interpret CELFIE count and reference evidence without accepting biological claims.

This module performs no I/O and never rewrites historical gates. Callers supply
the actual reference SHA256 and, when available, the retained compatibility
audit and its source identity. A passed count gate establishes count integrity;
it does not establish the original megakaryocyte fragment reference or accuracy.
"""
from copy import deepcopy
import math
import re


def _sha256(value):
    return isinstance(value, str) and re.fullmatch(r'[0-9a-f]{64}', value) is not None


def _nonnegative_integer(value):
    return type(value) is int and value >= 0


def _audit_reference_sha(evidence):
    """Read the reference role from the retained audit's explicit command."""
    script = evidence.get('script')
    argv = script.get('argv', []) if isinstance(script, dict) else []
    if not isinstance(argv, list) or argv.count('--reference') != 1:
        return None
    index = argv.index('--reference') + 1
    if index >= len(argv) or not isinstance(argv[index], str):
        return None
    inputs = evidence.get('inputs')
    record = inputs.get(argv[index]) if isinstance(inputs, dict) else None
    return record.get('sha256') if isinstance(record, dict) else None


def _contradiction(row):
    if not isinstance(row, dict):
        return False
    methylated = row.get('reconstructed_methylated')
    total = row.get('reconstructed_total')
    unmethylated = row.get('reconstructed_unmethylated')
    score = row.get('shipped_U_fraction')
    return (
        row.get('atlas_direction') == 'U'
        and row.get('marker_target') == 'Megakaryocytes'
        and _nonnegative_integer(row.get('atlas_row_1based'))
        and row['atlas_row_1based'] > 0
        and _nonnegative_integer(methylated)
        and _nonnegative_integer(total)
        and methylated == total > 0
        and type(unmethylated) is int and unmethylated == 0
        and type(score) in (int, float) and math.isfinite(score)
        and 0 < score <= 1
    )


def evaluate_reference_status(gate: dict, reference_sha: str, *, evidence=None,
                              evidence_source=None) -> dict:
    """Return conservative, separate count-integrity and reference conclusions.

    ``evidence`` supports the retained megakaryocyte compatibility-audit schema.
    Its reference and atlas hashes must match the supplied matrix and gate;
    otherwise the megakaryocyte identity remains unresolved. A contradiction
    refers to the same complete, strictly region-clipped observations, not to
    equivalence between methylation percentages and fragment fractions in general.

    Neither this audit nor the older count gates establish biological accuracy.
    Accordingly this policy has no automatic path to an accepted matched or
    biologically validated reference. Future positive evidence needs an explicit
    reviewed policy, rather than recognition of a free-text equivalence claim.
    """
    gate = gate if isinstance(gate, dict) else {}
    reasons = []
    valid_sha = _sha256(reference_sha)
    recorded_sha = gate.get('reference_counts_sha256')
    hashes_match = valid_sha and recorded_sha == reference_sha
    if gate.get('status') == 'passed' and hashes_match:
        integrity = 'passed'
    elif gate.get('status') not in (None, 'passed') or (
            valid_sha and _sha256(recorded_sha) and not hashes_match):
        integrity = 'failed'
        reasons.append('The count gate failed or does not identify the actual reference counts.')
    else:
        integrity = 'unverified'
        reasons.append('A passed count gate bound to the actual reference SHA256 is missing.')

    checked = gate.get('pat_cell_columns_checked')
    passed = gate.get('pat_cell_columns_passed')
    pat_verified = (integrity == 'passed' and checked == 39 and passed == 39
                    and type(checked) is int and type(passed) is int
                    and gate.get('unresolved_cells') == [])
    pat_status = {
        'status': 'verified_at_published_precision' if pat_verified else 'not_established',
        'checked': checked if _nonnegative_integer(checked) else None,
        'passed': passed if _nonnegative_integer(passed) else None,
        'scope': '39 PAT-backed columns; excludes megakaryocytes; numerical reconstruction only',
    }
    reference_status = 'unresolved'
    evidence_binding = 'not_supplied'
    conflicts = []
    if isinstance(evidence, dict):
        atlas_sha = gate.get('atlas_sha256')
        if not (valid_sha and _audit_reference_sha(evidence) == reference_sha
                and _sha256(atlas_sha) and evidence.get('fixed_atlas_sha256') == atlas_sha):
            evidence_binding = 'hash_mismatch_or_missing_identity'
            reasons.append('Compatibility evidence is not bound to this reference and atlas.')
        elif (evidence.get('status') != 'compatibility_failure'
              or evidence.get('MK_component_exactly_reproduces_all_979_reference_METH_DEPTH_pairs') is not True
              or evidence.get('MK_region_names_and_CpG_interval_lengths_match_component') is not True):
            evidence_binding = 'unsupported_or_incomplete_evidence'
            reasons.append('The supplied audit does not establish the retained incompatibility proof.')
        else:
            rows = evidence.get('incompatible_MK_target_regions')
            if isinstance(rows, list) and rows and all(_contradiction(row) for row in rows):
                evidence_binding = 'matched'
                reference_status = 'incompatible_with_shipped_uxm'
                conflicts = [row['atlas_row_1based'] for row in rows]
                reasons.append('The bound audit demonstrates incompatible megakaryocyte observations at the recorded regions.')
            else:
                evidence_binding = 'unsupported_or_incomplete_evidence'
                reasons.append('The supplied audit lacks internally consistent contradictory region records.')
    if reference_status == 'unresolved':
        reasons.append('Original megakaryocyte fragment-reference identity is not established.')
    reasons.append('Count integrity and numerical convergence do not establish biological accuracy.')
    return {
        'schema_version': 1,
        'reference_counts_sha256': reference_sha,
        'atlas_sha256': gate.get('atlas_sha256'),
        'count_integrity': integrity,
        'count_gate_hash_matches': hashes_match,
        'pat_backed_identity': pat_status,
        'reference_status': reference_status,
        'megakaryocyte_identity': reference_status,
        'incompatible_atlas_rows_1based': conflicts,
        'compatibility_evidence': {
            'binding': evidence_binding,
            'source': deepcopy(evidence_source),
            'scope': 'Same complete observations under strict regional clipping',
        },
        'matched_uxm_reference': False,
        'biological_validation': 'not_established',
        'diagnostic_only': True,
        'reasons': reasons,
    }
