#!/usr/bin/env python3
"""Count actual CpG observations and native UXM fragment classes in complete PAT streams."""
import argparse
import sys
import concurrent.futures
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading

COUNTER_CPP = r'''
#include <algorithm>
#include <charconv>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
struct Region { std::string chr; int64_t s,e,m=0,t=0; };
int64_t integer(const std::string& x) { int64_t n; auto r=std::from_chars(x.data(),x.data()+x.size(),n); if(r.ec!=std::errc() || r.ptr!=x.data()+x.size()) throw std::runtime_error("Invalid integer"); return n; }
std::vector<std::string> fields(const std::string& s) { std::vector<std::string> v; size_t a=0,b; while((b=s.find('\t',a))!=std::string::npos) { v.push_back(s.substr(a,b-a)); a=b+1; } v.push_back(s.substr(a)); return v; }
void add(int64_t& dst,int64_t mult,int64_t n) { if(n && mult>(std::numeric_limits<int64_t>::max()-dst)/n) throw std::runtime_error("Count overflow"); dst+=mult*n; }
int main(int argc,char**argv) { try {
 if(argc!=4) throw std::runtime_error("Usage: counter blocks mt.tsv stats.json");
 std::ifstream blocks(argv[1]); if(!blocks) throw std::runtime_error("Cannot open blocks");
 std::vector<Region> regions; std::string line;
 while(std::getline(blocks,line)) { auto v=fields(line); if(v.size()<5) throw std::runtime_error("Invalid block"); regions.push_back({v[0],integer(v[3]),integer(v[4])}); }
 if(regions.empty()) throw std::runtime_error("No blocks");
 for(size_t i=0;i<regions.size();++i) if(regions[i].s<1 || regions[i].e<=regions[i].s || (i&&regions[i].s<regions[i-1].e)) throw std::runtime_error("Overlapping or unsorted blocks");
 size_t ri=0; int64_t previous=0,rows=0,selected=0,maxlen=0,fragments=0;
 std::ios::sync_with_stdio(false); std::cin.tie(nullptr);
 while(std::getline(std::cin,line)) {
  if(line.empty() || line[0]=='#') continue;
  auto v=fields(line); if(v.size()<4) throw std::runtime_error("PAT row has fewer than four columns");
  int64_t s=integer(v[1]), n=integer(v[3]); const std::string& pat=v[2];
  if(s<1 || s<previous || n<1 || n>std::numeric_limits<int32_t>::max() || pat.empty()) throw std::runtime_error("Invalid or unsorted PAT");
  previous=s; ++rows; add(fragments,n,1); maxlen=std::max(maxlen,(int64_t)pat.size());
  for(char c:pat) if(c!='C'&&c!='T'&&c!='.') throw std::runtime_error("PAT contains unsupported character");
  const int64_t e=s+(int64_t)pat.size();
  while(ri<regions.size()&&regions[ri].e<=s) ++ri;
  bool chosen=false;
  for(size_t j=ri;j<regions.size()&&regions[j].s<e;++j) {
   auto& r=regions[j]; if(r.e<=s) continue;
   if(r.chr!=v[0]) throw std::runtime_error("PAT chromosome disagrees with global CpG index");
   int64_t nc=0,nt=0;
   for(int64_t k=std::max(s,r.s)-s;k<std::min(e,r.e)-s;++k) { nc+=(pat[k]=='C'); nt+=(pat[k]=='T'); }
   add(r.m,n,nc); add(r.t,n,nc+nt); chosen=true;
  }
  if(chosen) { std::cout<<line<<'\n'; ++selected; }
 }
 if(!std::cin.eof()) throw std::runtime_error("Input stream read failed");
 std::ofstream out(argv[2]); if(!out) throw std::runtime_error("Cannot write count output");
 for(const auto&r:regions) out<<r.m<<'\t'<<r.t<<'\n';
 out.close(); if(!out) throw std::runtime_error("Count output write failed");
 std::ofstream stats(argv[3]); stats<<"{\"pat_rows\":"<<rows<<",\"selected_rows\":"<<selected<<",\"represented_fragments\":"<<fragments<<",\"maximum_pattern_cpgs\":"<<maxlen<<"}\n";
 if(!stats) throw std::runtime_error("Stats output write failed");
 return 0;
 }catch(const std::exception&e) { std::cerr<<"PAT counter: "<<e.what()<<'\n'; return 1; }
}
'''


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    stage = path.with_name(path.name + '.partial')
    stage.write_text(json.dumps(value, indent=2) + '\n')
    stage.replace(path)


def compile_counter(output):
    """Build an immutable native overlap counter; return provenance and executable."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    source = output / 'pat_overlap_counter.cpp'
    expected = hashlib.sha256(COUNTER_CPP.encode()).hexdigest()
    binary = output / 'pat_overlap_counter'
    manifest = output / 'compiler.json'
    if manifest.exists():
        previous = json.loads(manifest.read_text())
        if previous['source_sha256'] != expected or sha(source) != expected or sha(binary) != previous['binary_sha256']:
            raise ValueError('Compiled counter identity differs from the existing run')
        return binary, previous
    if source.exists() and sha(source) != expected:
        raise ValueError('Partial counter source differs')
    source.write_text(COUNTER_CPP)
    compiler = shutil.which('g++')
    if not compiler:
        raise ValueError('g++ is required to build the small streaming counter')
    command = [compiler, '-O3', '-std=c++17', '-Wall', '-Wextra', str(source), '-o', str(binary)]
    started = utc()
    subprocess.run(command, check=True)
    provenance = {'started_utc': started, 'completed_utc': utc(), 'compiler': compiler,
                  'compiler_version': subprocess.check_output([compiler, '--version'], text=True).splitlines()[0],
                  'argv': command, 'source_sha256': sha(source), 'binary_sha256': sha(binary)}
    atomic_json(manifest, provenance)
    return binary, provenance


def atlas_blocks(atlas, output):
    with Path(atlas).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    content = ''.join('\t'.join(row[key] for key in ('chr', 'start', 'end', 'startCpG', 'endCpG')) + '\n' for row in rows)
    output = Path(output)
    if output.exists() and output.read_text() != content:
        raise ValueError('Existing count blocks differ from atlas')
    output.write_text(content)
    return rows


def stat_identity(path):
    stat = Path(path).stat()
    return {'path': str(Path(path).resolve()), 'size_bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def count_pat(pat, atlas, output, tools_dir, *, compiled=None):
    """Read a complete gzip PAT once, hash its compressed bytes and count exact observations.

    The native counter clips each pattern to [startCpG,endCpG), counts every C/T
    times PAT multiplicity, and forwards overlapping rows unchanged to pinned
    wgbstools homog with l4. Dots contribute neither methylated nor total counts.
    """
    pat, atlas, output, tools_dir = map(Path, (pat, atlas, output, tools_dir))
    output.mkdir(parents=True, exist_ok=True)
    native = tools_dir / 'share/wgbs_tools/src/homog/homog'
    binary, compiler = compiled or compile_counter(output / 'counter')
    blocks = output / 'blocks.tsv'
    rows = atlas_blocks(atlas, blocks)
    identity = stat_identity(pat)
    settings = {'atlas_sha256': sha(atlas), 'helper_sha256': sha(__file__),
                'counter_source_sha256': compiler['source_sha256'], 'counter_binary_sha256': compiler['binary_sha256'],
                'native_homog_sha256': sha(native), 'minimum_fragment_cpgs': 4,
                'fragment_ranges': '0,0.251,0.75,1', 'fragment_inclusive': False,
                'cpg_counting': 'all C/T observations, no fragment minimum, each multiplied by PAT count',
                'input': identity}
    execution = output / 'execution.json'
    result_path = output / 'counts.tsv'
    prior = json.loads(execution.read_text()) if execution.exists() else None
    if prior and prior['settings'] != settings:
        raise ValueError(f'Cannot resume changed count inputs or settings: {output}')
    if prior and prior.get('status') == 'completed':
        if sha(result_path) != prior['counts_sha256']:
            raise ValueError(f'Cached count checksum differs: {output}')
        return prior
    provenance = prior or {'started_utc': utc(), 'settings': settings, 'attempts': []}
    provenance['status'] = 'running'
    provenance['attempts'].append({'started_utc': utc()})
    atomic_json(execution, provenance)
    mt = output / 'mt.partial.tsv'
    stats = output / 'stats.partial.json'
    native_output = output / 'uxm.partial.tsv'
    compressed_hash = hashlib.sha256()
    feed_errors = []
    commands = [['gzip', '-dc'], [str(binary), str(blocks), str(mt), str(stats)],
                [str(native), '-b', str(blocks), '-l', '4', '-r', '0,0.251,0.75,1', '-n', pat.name]]
    provenance['attempts'][-1]['commands'] = commands
    atomic_json(execution, provenance)
    try:
        with (output / 'counter.log').open('ab') as log, native_output.open('wb') as native_handle:
            decompressor = subprocess.Popen(commands[0], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log)
            counter = subprocess.Popen(commands[1], stdin=decompressor.stdout, stdout=subprocess.PIPE, stderr=log)
            decompressor.stdout.close()
            homog = subprocess.Popen(commands[2], stdin=counter.stdout, stdout=native_handle, stderr=log)
            counter.stdout.close()
            def feed():
                try:
                    with pat.open('rb') as handle:
                        for data in iter(lambda: handle.read(1024 * 1024), b''):
                            compressed_hash.update(data)
                            decompressor.stdin.write(data)
                except BaseException as error:
                    feed_errors.append(repr(error))
                finally:
                    try:
                        decompressor.stdin.close()
                    except BrokenPipeError:
                        pass
            feeder = threading.Thread(target=feed, daemon=True)
            feeder.start()
            codes = [homog.wait(), counter.wait(), decompressor.wait()]
            feeder.join()
        if feed_errors or any(code != 0 for code in codes):
            raise ValueError(f'PAT counting pipeline failed: codes={codes}; feeder={feed_errors}; see {output / "counter.log"}')
        if stat_identity(pat) != identity:
            raise ValueError('PAT changed while counting')
        mt_values = [list(map(int, line.split())) for line in mt.read_text().splitlines()]
        ux_values = [list(map(int, line.split())) for line in native_output.read_text().splitlines()]
        if len(mt_values) != len(rows) or len(ux_values) != len(rows):
            raise ValueError('Native count row number differs from atlas')
        stage = output / 'counts.partial.tsv'
        with stage.open('w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['chr', 'start', 'end', 'startCpG', 'endCpG', 'methylated', 'total', 'U', 'X', 'M'])
            for row, pair, triplet in zip(rows, mt_values, ux_values):
                if len(pair) != 2 or len(triplet) != 3 or any(value < 0 for value in pair + triplet) or pair[0] > pair[1]:
                    raise ValueError('Invalid count values')
                writer.writerow([row[key] for key in ('chr', 'start', 'end', 'startCpG', 'endCpG')] + pair + triplet)
        stage.replace(result_path)
        provenance.update(status='completed', completed_utc=utc(), counts_sha256=sha(result_path),
                          compressed_pat_sha256=compressed_hash.hexdigest(), stream=json.loads(stats.read_text()))
        provenance['attempts'][-1]['completed_utc'] = utc()
        atomic_json(execution, provenance)
        return provenance
    except BaseException as error:
        provenance['status'] = 'failed'
        provenance['attempts'][-1].update(failed_utc=utc(), error=repr(error))
        atomic_json(execution, provenance)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pat', type=Path, required=True)
    parser.add_argument('--atlas', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tools-dir', type=Path, default=Path('/path/to/hpc-home/.conda/envs/nextflow'))
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from common import output_guard, PROJECT
    revision = output_guard.require_clean_source(PROJECT)
    dt.datetime.strptime(args.output.name[:16], '%Y-%m-%d_%H-%M')
    if len(args.output.name) < 18 or args.output.name[16] != '-':
        parser.error('Output name must start YYYY-MM-DD_HH-MM-')
    output_guard.reserve_path(args.output, {'project_revision': revision,
        'source_sha256': sha(__file__), 'pat': stat_identity(args.pat),
        'atlas_sha256': sha(args.atlas), 'tools_dir': str(args.tools_dir.resolve())})
    output_guard.update(args.output, status='running')
    try:
        result = count_pat(args.pat, args.atlas, args.output, args.tools_dir)
        output_guard.finish(args.output, True)
        print(json.dumps(result, indent=2))
    except BaseException:
        output_guard.finish(args.output, False)
        raise


if __name__ == '__main__':
    main()
