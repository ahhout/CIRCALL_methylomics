"""Regression checks separating count integrity from reference interpretation."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.celfie_reference_status import evaluate_reference_status


class ReferenceStatus(unittest.TestCase):
    def setUp(self):
        self.reference_sha = 'a' * 64
        self.atlas_sha = 'b' * 64
        self.gate = {
            'status': 'passed', 'reference_counts_sha256': self.reference_sha,
            'atlas_sha256': self.atlas_sha, 'pat_cell_columns_checked': 39,
            'pat_cell_columns_passed': 39, 'unresolved_cells': [],
        }
        self.evidence = {
            'status': 'compatibility_failure', 'fixed_atlas_sha256': self.atlas_sha,
            'script': {'argv': ['audit.py', '--reference', '/historic/reference.tsv']},
            'inputs': {'/historic/reference.tsv': {'sha256': self.reference_sha}},
            'MK_component_exactly_reproduces_all_979_reference_METH_DEPTH_pairs': True,
            'MK_region_names_and_CpG_interval_lengths_match_component': True,
            'incompatible_MK_target_regions': [{
                'atlas_row_1based': 96, 'atlas_direction': 'U',
                'marker_target': 'Megakaryocytes', 'shipped_U_fraction': 0.903,
                'reconstructed_methylated': 24, 'reconstructed_total': 24,
                'reconstructed_unmethylated': 0,
            }],
        }

    def test_passed_39_column_gate_never_implies_matched_40_cell_reference(self):
        result = evaluate_reference_status(self.gate, self.reference_sha)
        self.assertEqual(result['count_integrity'], 'passed')
        self.assertEqual(result['pat_backed_identity']['status'], 'verified_at_published_precision')
        self.assertEqual(result['reference_status'], 'unresolved')
        self.assertFalse(result['matched_uxm_reference'])
        self.assertTrue(result['diagnostic_only'])
        self.assertEqual(result['biological_validation'], 'not_established')

    def test_matching_audit_identifies_incompatibility_and_retains_source(self):
        source = {'path': '/audit/evidence.json', 'sha256': 'c' * 64}
        result = evaluate_reference_status(self.gate, self.reference_sha,
                                           evidence=self.evidence, evidence_source=source)
        self.assertEqual(result['reference_status'], 'incompatible_with_shipped_uxm')
        self.assertEqual(result['incompatible_atlas_rows_1based'], [96])
        self.assertEqual(result['compatibility_evidence']['binding'], 'matched')
        self.assertEqual(result['compatibility_evidence']['source'], source)
        self.assertTrue(result['diagnostic_only'])
        self.assertFalse(result['matched_uxm_reference'])
        source['path'] = '/changed'
        self.assertEqual(result['compatibility_evidence']['source']['path'], '/audit/evidence.json')

    def test_changed_reference_is_not_covered_by_old_gate_or_old_audit(self):
        result = evaluate_reference_status(self.gate, 'd' * 64, evidence=self.evidence)
        self.assertEqual(result['count_integrity'], 'failed')
        self.assertEqual(result['reference_status'], 'unresolved')
        self.assertEqual(result['pat_backed_identity']['status'], 'not_established')
        self.assertTrue(result['diagnostic_only'])

    def test_wrong_atlas_or_reference_evidence_fails_closed(self):
        for target in ('reference', 'atlas'):
            with self.subTest(target=target):
                evidence = deepcopy(self.evidence)
                if target == 'reference':
                    evidence['inputs']['/historic/reference.tsv']['sha256'] = 'd' * 64
                else:
                    evidence['fixed_atlas_sha256'] = 'd' * 64
                result = evaluate_reference_status(self.gate, self.reference_sha, evidence=evidence)
                self.assertEqual(result['count_integrity'], 'passed')
                self.assertEqual(result['reference_status'], 'unresolved')
                self.assertEqual(result['compatibility_evidence']['binding'], 'hash_mismatch_or_missing_identity')

    def test_incomplete_or_malformed_proof_does_not_assert_incompatibility(self):
        for field, value in [('reconstructed_methylated', 23), ('reconstructed_total', 0),
                             ('shipped_U_fraction', 0), ('shipped_U_fraction', float('nan')),
                             ('shipped_U_fraction', True), ('reconstructed_unmethylated', False),
                             ('atlas_direction', 'M')]:
            with self.subTest(field=field, value=value):
                evidence = deepcopy(self.evidence)
                evidence['incompatible_MK_target_regions'][0][field] = value
                result = evaluate_reference_status(self.gate, self.reference_sha, evidence=evidence)
                self.assertEqual(result['reference_status'], 'unresolved')
                self.assertTrue(result['diagnostic_only'])

    def test_free_text_or_new_claims_cannot_promote_unestablished_reference(self):
        self.gate.update(megakaryocyte_original_uxm_pattern_equivalence='verified',
                         matched_uxm_reference=True, biological_validation='passed',
                         diagnostic_only=False)
        result = evaluate_reference_status(self.gate, self.reference_sha)
        self.assertEqual(result['reference_status'], 'unresolved')
        self.assertFalse(result['matched_uxm_reference'])
        self.assertEqual(result['biological_validation'], 'not_established')
        self.assertTrue(result['diagnostic_only'])

    def test_combined_gate_cannot_inherit_pat_identity_from_description(self):
        self.gate.pop('pat_cell_columns_checked')
        self.gate.pop('pat_cell_columns_passed')
        self.gate.update(reference_scope='39 numerically checked PAT pools',
                         reference_counts={'path': '/new/reference.tsv', 'sha256': self.reference_sha})
        result = evaluate_reference_status(self.gate, self.reference_sha, evidence=self.evidence)
        self.assertEqual(result['count_integrity'], 'passed')
        self.assertEqual(result['pat_backed_identity']['status'], 'not_established')
        self.assertEqual(result['reference_status'], 'incompatible_with_shipped_uxm')

    def test_missing_gate_hash_and_failed_gate_never_verify_counts(self):
        result = evaluate_reference_status({}, '')
        self.assertEqual(result['count_integrity'], 'unverified')
        self.assertFalse(result['count_gate_hash_matches'])
        self.gate['status'] = 'failed'
        result = evaluate_reference_status(self.gate, self.reference_sha, evidence=self.evidence)
        self.assertEqual(result['count_integrity'], 'failed')
        self.assertEqual(result['pat_backed_identity']['status'], 'not_established')

    def test_audit_must_identify_reference_role_unambiguously(self):
        for argv in ([], ['audit.py', '--reference'],
                     ['audit.py', '--reference', '/historic/reference.tsv', '--reference', '/other']):
            with self.subTest(argv=argv):
                evidence = deepcopy(self.evidence)
                evidence['script']['argv'] = argv
                result = evaluate_reference_status(self.gate, self.reference_sha, evidence=evidence)
                self.assertEqual(result['reference_status'], 'unresolved')

    def test_retained_historical_evidence_schema_is_supported_without_mutation(self):
        project = Path(__file__).resolve().parents[3]
        evidence_path = project / ('workflows/deconflow/runs/'
            '2026-10-26_13-27-uxm-celfie-matched-reference-v1/megakaryocyte_reference_compatibility.json')
        original = evidence_path.read_bytes()
        evidence = json.loads(original)
        gate = deepcopy(self.gate)
        argv = evidence['script']['argv']
        reference_path = argv[argv.index('--reference') + 1]
        gate['reference_counts_sha256'] = evidence['inputs'][reference_path]['sha256']
        gate['atlas_sha256'] = evidence['fixed_atlas_sha256']
        result = evaluate_reference_status(gate, gate['reference_counts_sha256'], evidence=evidence,
                                           evidence_source={'path': str(evidence_path)})
        self.assertEqual(result['reference_status'], 'incompatible_with_shipped_uxm')
        self.assertEqual(result['incompatible_atlas_rows_1based'], [96, 275])
        self.assertEqual(evidence_path.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
