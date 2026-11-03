#!/usr/bin/env python3
"""Bind the completed reference and patient count preparations for native CELFIE."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import ATLAS_SHA, WORKFLOW, atlas_rows, identity, read_json, rows, sha, utc, write_json
from celfie_patients import count_table


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-preparation', type=Path, required=True)
    p.add_argument('--patient-preparation', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from common import output_guard, PROJECT
    output_guard.require_clean_source(PROJECT)
    output_guard.canonical_parent(a.output.absolute().parent)
    started = utc()
    if a.output.exists():
        p.error('Combined gate exists; preserve it rather than overwriting')
    refgate_path = a.reference_preparation/'validation_gate.json'
    patgate_path = a.patient_preparation/'validation_gate.json'
    refgate, patgate = read_json(refgate_path), read_json(patgate_path)
    for label, gate in [('reference', refgate), ('patient', patgate)]:
        if gate.get('status') != 'passed' or gate.get('atlas_sha256') != ATLAS_SHA:
            raise ValueError(f'The {label} count preparation has not passed for the fixed atlas')
    reference, patients = a.reference_preparation/'reference.counts.tsv', a.patient_preparation/'patient.counts.tsv'
    if sha(reference) != refgate['reference_counts_sha256'] or sha(patients) != patgate['patient_counts_sha256']:
        raise ValueError('Count matrix differs from its preparation gate')
    atlas, cells = atlas_rows()
    coordinates = [(r['chr'],r['start'],r['end']) for r in atlas]
    count_table(reference, coordinates, cells)
    _, samples = count_table(patients, coordinates)
    expected = [r['sample'] for r in rows(WORKFLOW/'config/samples.pats.tsv')]
    if samples != expected or patgate['sample_order'] != expected:
        raise ValueError('Patient counts must retain all twelve approved IDs in recorded order')
    if refgate.get('pat_cell_columns_checked') != 39 or refgate.get('pat_cell_columns_passed') != 39:
        raise ValueError('All 39 PAT-based reference cell columns must pass the pooling check')
    document = dict(started_utc=started, completed_utc=utc(), status='passed',
        atlas_sha256=ATLAS_SHA, regions=979, cells=cells, samples=samples,
        reference_counts_sha256=sha(reference), patient_counts_sha256=sha(patients),
        reference_counts=identity(reference, True), patient_counts=identity(patients, True),
        source_gates={'reference':identity(refgate_path,True), 'patients':identity(patgate_path,True)},
        helper=identity(__file__,True),
        reference_scope='39 numerically checked PAT pools plus documented BLUEPRINT megakaryocyte count reconstruction',
        megakaryocyte_fragment_equivalence='Not established; per-CpG count reconstruction does not recover original UXM fragment patterns',
        count_policy='Actual integer CpG observations; preserve all 979 regions and real zero depth; native CELFIE pseudocounts unchanged')
    # Exclusive creation prevents a second caller replacing an accepted gate.
    import json
    with a.output.open('x') as f:
        json.dump(document,f,indent=2,allow_nan=False)
        f.write('\n')
    print(a.output)


if __name__=='__main__':
    main()
