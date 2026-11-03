#!/usr/bin/env python3
"""Freeze a bounded paired-read prefix from one production samplesheet entry."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / 'lib'))
import output_guard as guard

MAX_PAIRS = 10000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    require(path.stat().st_size <= 64 * 1024 * 1024, 'Smoke artifact exceeds bounded size')
    return hashlib.sha256(path.read_bytes()).hexdigest()


def production_row(project):
    with (project / 'samplesheets/samplesheet.csv').open(newline='') as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ['sample', 'fastq_1', 'fastq_2', 'genome'], 'Unexpected samplesheet columns')
        row = next(reader, None)
    require(row and row['sample'] and row['fastq_1'] and row['fastq_2'], 'First entry must contain paired FASTQs')
    return row


def source_state(row):
    state = {}
    for key in ('fastq_1', 'fastq_2'):
        p = Path(row[key]).resolve(strict=True)
        require(p.is_file(), 'FASTQ source must be a file')
        st = p.stat()
        state[key] = dict(path=str(p), size=st.st_size, mtime_ns=st.st_mtime_ns, inode=st.st_ino)
    return state


def record(handle):
    lines = [handle.readline() for _ in range(4)]
    require(all(lines) and lines[0].startswith(b'@') and lines[2].startswith(b'+'), 'Truncated or invalid FASTQ record')
    require(len(lines[1].rstrip()) == len(lines[3].rstrip()), 'FASTQ sequence/quality length mismatch')
    require(len(lines[1]) <= 2000, 'Unexpected long read in bounded short-read smoke input')
    return lines


def pair_id(lines):
    return lines[0].split()[0].removesuffix(b'/1').removesuffix(b'/2')


def prepare(project, name, pairs=MAX_PAIRS):
    project = project.resolve(strict=True)
    require(type(pairs) is int and 1 <= pairs <= MAX_PAIRS, 'Smoke requires 1–10000 read pairs')
    revision = guard.require_clean_source(project)
    row = production_row(project)
    sources = source_state(row)
    parent = project / 'runs/methylseq/smoke_inputs'
    require(parent.resolve() == parent, 'Smoke input parent cannot be an alias')
    parent.mkdir(parents=True, exist_ok=True)
    output = guard.reserve(parent, name, dict(project_commit=revision, pairs=pairs, sources=sources))
    guard.update(output, status='running')
    try:
        targets = [output / 'reads_R1.fastq.gz', output / 'reads_R2.fastq.gz']
        with gzip.open(sources['fastq_1']['path'], 'rb') as r1, gzip.open(sources['fastq_2']['path'], 'rb') as r2:
            with targets[0].open('xb') as f1, targets[1].open('xb') as f2:
                with gzip.GzipFile(filename='', mode='wb', fileobj=f1, mtime=0) as w1, gzip.GzipFile(filename='', mode='wb', fileobj=f2, mtime=0) as w2:
                    for _ in range(pairs):
                        a, b = record(r1), record(r2)
                        require(pair_id(a) == pair_id(b), 'FASTQ mates are out of sync')
                        w1.writelines(a); w2.writelines(b)
        require(source_state(row) == sources, 'Source FASTQ changed during extraction')
        subset = dict(row, fastq_1=str(targets[0]), fastq_2=str(targets[1]))
        sheet = output / 'samplesheet.csv'
        with sheet.open('x', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['sample', 'fastq_1', 'fastq_2', 'genome'])
            writer.writeheader(); writer.writerow(subset)
        manifest = dict(schema=1, selection='first production entry; first N paired records',
                        pairs=pairs, samplesheet=str(sheet), source_row=row, sources=sources,
                        production_samplesheet_sha256=sha(project / 'samplesheets/samplesheet.csv'),
                        artifact_sha256={p.name: sha(p) for p in [*targets, sheet]},
                        project_commit=revision, created_utc=guard.utc_now(),
                        scientific_settings='Production reference and parameters unchanged; scheduling bounded separately')
        destination = output / 'manifest.json'
        with destination.open('x') as handle:
            json.dump(manifest, handle, indent=2); handle.write('\n')
        guard.finish(output, True, manifest_sha256=sha(destination))
        return destination
    except BaseException as exc:
        guard.finish(output, False, failure=str(exc))
        raise


def validate(project, manifest_path):
    require(manifest_path.is_absolute() and ',' not in str(manifest_path)
            and manifest_path.resolve(strict=True) == manifest_path, 'Smoke manifest must be canonical')
    parent = project / 'runs/methylseq/smoke_inputs'
    require(manifest_path.parent.parent == parent and manifest_path.name == 'manifest.json', 'Smoke manifest outside project runtime')
    owner = guard.read_record(manifest_path.parent)
    require(owner['status'] == 'completed' and owner.get('manifest_sha256') == sha(manifest_path), 'Smoke input preparation incomplete or altered')
    manifest = json.loads(manifest_path.read_text())
    require(manifest.get('schema') == 1 and type(manifest.get('pairs')) is int and 1 <= manifest['pairs'] <= MAX_PAIRS, 'Invalid bounded smoke manifest')
    row = production_row(project)
    require(manifest['source_row'] == row and manifest['sources'] == source_state(row), 'Original smoke sources changed')
    require(manifest['production_samplesheet_sha256'] == sha(project / 'samplesheets/samplesheet.csv'), 'Production samplesheet changed')
    expected = {'reads_R1.fastq.gz', 'reads_R2.fastq.gz', 'samplesheet.csv'}
    require(set(manifest['artifact_sha256']) == expected, 'Smoke artifact inventory differs')
    for name, digest in manifest['artifact_sha256'].items():
        p = manifest_path.parent / name
        require(not p.is_symlink() and p.is_file() and sha(p) == digest, 'Smoke artifact changed')
    sheet = manifest_path.parent / 'samplesheet.csv'
    require(manifest['samplesheet'] == str(sheet), 'Smoke samplesheet escaped its input bundle')
    with sheet.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    require(rows == [dict(row, fastq_1=str(sheet.parent/'reads_R1.fastq.gz'), fastq_2=str(sheet.parent/'reads_R2.fastq.gz'))], 'Smoke sample identity or paired inputs differ')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-name', required=True)
    parser.add_argument('--pairs', type=int, default=MAX_PAIRS)
    args = parser.parse_args()
    print(prepare(Path(__file__).resolve().parents[1], args.run_name, args.pairs))
