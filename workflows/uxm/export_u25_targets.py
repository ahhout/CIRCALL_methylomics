#!/usr/bin/env python3
"""Export the pinned, complete hg38 UXM U25/l4 atlas for custom-panel review."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

WORKFLOW = Path(__file__).resolve().parent
PROJECT = WORKFLOW.parents[1]
sys.path.insert(0, str(WORKFLOW.parent))
from output_contract import output_guard
ATLAS_SHA256 = 'e58fbec4fe6200a311540c6e1a16be901a2d1e9ed9d74b6a57888f2feb635f51'
DICTIONARY_SHA256 = '28180d18667dc0a6e66dd4250593e03a8bf10be50e9212aefd4e9282afea8590'
UXM_COMMIT = '8d0bb456742ece21802bc8e651e080869c67c944'
EXPECTED = {'intervals':979, 'target_bp':312012, 'unique_cpgs':8861, 'cell_types':40, 'excluded_intervals':0}


def utc():
    return datetime.now(timezone.utc).isoformat()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def dump_json(path, value):
    with path.open('w') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def write_table(path, fields, rows):
    with path.open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def chrom_key(chrom):
    suffix = chrom.removeprefix('chr')
    require(suffix.isdigit() or suffix in ('X','Y'), f'Unexpected atlas chromosome: {chrom}')
    return int(suffix) if suffix.isdigit() else {'X':23, 'Y':24}[suffix]


def execute(args, started):
    stamp = datetime.fromisoformat(started).strftime('%Y-%m-%d_%H-%M')
    datetime.strptime(stamp, '%Y-%m-%d_%H-%M')
    out = args.output_root.absolute()/f'{stamp}-uxm-u25-targets-v1'
    require(not out.exists(), f'Refusing existing output directory: {out}')
    require(not args.output_root.exists() or args.output_root.is_dir(), 'Output root is not a directory')
    atlas, dictionary, tabix = args.atlas.resolve(), args.dictionary.resolve(), args.tabix.resolve()
    reverse = dictionary.parent/'rev.CpG.bed.gz'
    reverse_index = Path(str(reverse)+'.tbi')
    for path in (atlas,dictionary,tabix,reverse,reverse_index):
        require(path.is_file(), f'Missing required input: {path}')
    require(reverse.resolve() == dictionary, 'Reverse CpG index does not refer to the same dictionary')
    require(sha(atlas) == ATLAS_SHA256, 'Atlas SHA-256 differs from the pinned production U25 atlas')
    require(sha(dictionary) == DICTIONARY_SHA256, 'hg38 CpG dictionary SHA-256 mismatch')
    with atlas.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        rows = list(reader)
    require(fields[:8] == ['chr','start','end','startCpG','endCpG','target','name','direction'],
            'Unexpected atlas metadata columns')
    cell_types = fields[8:]
    require(len(rows) == 979 and len(cell_types) == 40, 'Expected 979 markers and 40 reference cell types')
    require(len({r['name'] for r in rows}) == 979, 'Duplicate source marker names')
    require({r['direction'] for r in rows} == {'U'}, 'Unexpected atlas marker direction')
    require(all(r['target'] in cell_types for r in rows), 'Unknown target cell-type label')
    for row in rows:
        require(int(row['endCpG'])-int(row['startCpG']) >= 4,
                f"Published marker cannot support l4: {row['name']}; refusing unapproved exclusions")
        for cell in cell_types:
            if row[cell] == 'NA':
                continue
            value = float(row[cell])
            require(math.isfinite(value) and 0 <= value <= 1, f'Invalid U fraction: {row["name"]}/{cell}')
    # This index queries the third (CpG-index) column, rather than genomic bases.
    queries = [f"{r['chr']}:{r['startCpG']}-{int(r['endCpG'])-1}" for r in rows]
    result = subprocess.run([str(tabix),str(reverse),*queries],check=True,capture_output=True,text=True)
    positions = {}
    for line in result.stdout.splitlines():
        chrom, pos, index = line.split()
        index, pos = int(index), int(pos)
        require(index not in positions or positions[index] == (chrom,pos), 'Conflicting reference CpG mappings')
        positions[index] = (chrom,pos)
    for row in rows:
        first, stop = int(row['startCpG']), int(row['endCpG'])
        require(all(i in positions and positions[i][0] == row['chr'] for i in range(first,stop)),
                f"Missing reference CpGs: {row['name']}")
        # Every published source start equals the 1-based first cytosine; source end is the final G.
        require(int(row['start']) == positions[first][1] and int(row['end']) == positions[stop-1][1]+1,
                f"Unverified genomic coordinate convention: {row['name']}")
    annotated = []
    indexed = sorted(enumerate(rows,1),key=lambda item:(chrom_key(item[1]['chr']),int(item[1]['start']),int(item[1]['end'])))
    for number,(source_row,row) in enumerate(indexed,1):
        first,stop = int(row['startCpG']),int(row['endCpG'])
        start,end = int(row['start'])-1,int(row['end'])
        require(0 <= start < end, 'Invalid converted target interval')
        marker = dict(chrom=row['chr'],start_0based=start,end_exclusive=end,
                      target_id=f'UXM_U25_{number:04d}',length_bp=end-start,cpg_count=stop-first,
                      cpg_indices=';'.join(map(str,range(first,stop))),
                      cytosine_positions_1based=';'.join(str(positions[i][1]) for i in range(first,stop)),
                      atlas_row_1based=source_row,atlas_marker_name=row['name'],
                      atlas_target_cell_type=row['target'],atlas_direction=row['direction'],
                      atlas_start_1based=int(row['start']),atlas_end_inclusive=int(row['end']),
                      atlas_start_cpg_index=first,atlas_end_cpg_index_exclusive=stop,
                      missing_reference_cell_types=';'.join(cell for cell in cell_types if row[cell]=='NA') or '.',
                      missing_reference_count=sum(row[cell]=='NA' for cell in cell_types),
                      minimum_fragment_cpgs=4,max_methylated_fraction_for_U='0.25')
        marker.update({f'atlas_U_fraction__{cell}':row[cell] for cell in cell_types})
        annotated.append(marker)
    require(len({(r['chrom'],r['start_0based'],r['end_exclusive']) for r in annotated}) == 979,
            'Duplicate converted target intervals')
    for left,right in zip(annotated,annotated[1:]):
        require(left['chrom'] != right['chrom'] or left['end_exclusive'] <= right['start_0based'],
                'Overlapping published marker intervals; refusing silent merging')
    selected_indices = [i for row in rows for i in range(int(row['startCpG']),int(row['endCpG']))]
    totals = dict(intervals=len(annotated),target_bp=sum(r['length_bp'] for r in annotated),
                  unique_cpgs=len(set(selected_indices)),cell_types=len(cell_types),excluded_intervals=0)
    require(totals == EXPECTED and len(selected_indices) == 8861, f'Unexpected export totals: {totals}')
    input_hashes = {str(path):sha(path) for path in (atlas,reverse_index,tabix,Path(__file__).resolve())}
    input_hashes[str(dictionary)] = DICTIONARY_SHA256
    revision = output_guard.require_clean_source(PROJECT)
    args.output_root.mkdir(parents=True,exist_ok=True)
    output_guard.reserve_path(out, {'workflow': str(Path(__file__)), 'input_sha256': input_hashes, 'project_revision': revision})
    provenance = dict(started_utc=started,status='running',argv=sys.argv,genome='hg38',candidate_export=True,
                      output_dir=str(out),input_sha256=input_hashes,observed_totals=totals,
                      source_missing_profiles={'markers_with_any_NA':sum(any(r[c]=='NA' for c in cell_types) for r in rows),
                                               'markers_with_all_NA':sum(all(r[c]=='NA' for c in cell_types) for r in rows),
                                               'NA_values':sum(r[c]=='NA' for r in rows for c in cell_types),
                                               'handling':'Preserved literally as published and flagged; no imputation or exclusion'},
                      source_atlas=str(atlas),source_commit=UXM_COMMIT,
                      source_url=f'https://raw.githubusercontent.com/nloyfer/UXM_deconv/{UXM_COMMIT}/supplemental/Atlas.U25.l4.hg38.full.tsv',
                      coordinate_validation={'markers_verified':979,'source_start':'1-based first CpG cytosine',
                                             'source_end':'1-based inclusive final CpG guanine',
                                             'conversion':'BED start = source start - 1; BED end = source end',
                                             'method':'Verify every source CpG index and both genomic endpoints against the pinned hg38 dictionary'},
                      conventions={'targets':'0-based half-open; no header in BED',
                                   'cpg_indices':'wgbstools hg38 indices; source endCpG is exclusive',
                                   'cytosine_positions':'1-based top-strand cytosines',
                                   'atlas_U_fraction_columns':'Original unmethylated-fragment proportions, copied without numerical rewriting',
                                   'selection':'All published markers and all 40 cell types; chrX/chrY retained as published; no merging or trimming',
                                   'fragment_definition':'Official l4 atlas: at least four informative CpGs; U means at most 25% methylated CpGs',
                                   'exclusions':'None; empty excluded_intervals.tsv records that no PCDH ten-CpG filter was applied',
                                   'limitations':'Candidate targets; no probe padding, tiling, sequence design or assay validation'},
                      software={'python':sys.version,'dependencies':'Python standard library and existing tabix'})
    dump_json(out/'provenance.json',provenance)
    try:
        with (out/'targets.hg38.bed').open('x') as handle:
            for row in annotated:
                handle.write(f"{row['chrom']}\t{row['start_0based']}\t{row['end_exclusive']}\t{row['target_id']}\n")
        write_table(out/'targets.annotated.tsv',list(annotated[0]),annotated)
        write_table(out/'excluded_intervals.tsv',
                    ['atlas_row_1based','atlas_marker_name','chrom','source_start_1based','source_end_inclusive','reason'],[])
        shutil.copyfile(__file__,out/Path(__file__).name)
        with (out/'targets.annotated.tsv').open() as handle:
            written = list(csv.DictReader(handle,delimiter='\t'))
        bed_rows = [line.split('\t') for line in (out/'targets.hg38.bed').read_text().splitlines()]
        require(bed_rows == [[r['chrom'],r['start_0based'],r['end_exclusive'],r['target_id']] for r in written],
                'Written BED and annotated table disagree')
        require(len(written) == 979 and len((out/'excluded_intervals.tsv').read_text().splitlines()) == 1,
                'Incomplete output or unintended exclusions')
        provenance.update(status='completed',completed_utc=utc(),
                          checks_passed=['Pinned atlas and reference dictionary SHA-256 hashes match',
                                         'Every marker CpG index and endpoint validated against hg38',
                                         '979 unique non-overlapping regions retained; no exclusions',
                                         'All 40 reference cell-type profiles retained',
                                         '312012 bp and 8861 unique CpGs reconcile',
                                         'Written BED and annotated TSV coordinates/IDs agree'],
                          output_sha256={p.name:sha(p) for p in sorted(out.iterdir()) if p.name!='provenance.json'})
        dump_json(out/'provenance.json',provenance)
        output_guard.finish(out, provenance['status'] == 'completed')
    except Exception as error:
        provenance.update(status='failed',failed_utc=utc(),error=str(error))
        dump_json(out/'provenance.json',provenance)
        output_guard.finish(out, provenance['status'] == 'completed')
        raise
    print(out)
    print(json.dumps(totals,indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--atlas',type=Path,default=WORKFLOW/'reference/Atlas.U25.l4.hg38.full.tsv')
    parser.add_argument('--dictionary',type=Path,default=PROJECT/'references/wgbstools_hg38/CpG.bed.gz')
    parser.add_argument('--tabix',type=Path,default=Path(sys.executable).parent/'tabix')
    parser.add_argument('--output-root',type=Path,default=WORKFLOW/'pannel_cpg_list')
    args = parser.parse_args()
    try:
        execute(args,utc())
    except (ValueError,OSError,KeyError,TypeError,subprocess.CalledProcessError) as error:
        parser.exit(1,f'UXM U25 target export failed: {error}\n')


if __name__ == '__main__':
    main()
