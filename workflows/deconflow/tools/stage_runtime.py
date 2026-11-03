#!/usr/bin/env python3
"""Stage immutable public runtime resources; only this step accesses the network."""
import argparse
import sys
import os
import shutil
import tarfile
import urllib.request
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import *


def cached_download(path, url, expected):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        # Exclusive creation preserves any previous download or failure evidence.
        with urllib.request.urlopen(url, timeout=120) as response, path.open('xb') as out:
            shutil.copyfileobj(response, out)
    if sha(path) != expected:
        raise ValueError(f'Download checksum mismatch: {path}; retained for inspection')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    root = new_run(args.output)
    try:
        (root / 'images').mkdir()
        (root / 'tmp').mkdir()
        shutil.copytree(WORKFLOW / 'config', root / 'config')
        shutil.copy2(WORKFLOW / 'environment.yml', root)
        source = read_json(root / 'config/source.lock.json')
        archive = WORKFLOW / 'software/upstream.tar.gz'
        cached_download(archive, source['archive_url'], source['archive_sha256'])
        cached_download(WORKFLOW / 'software/loyfer_supplement.xlsx',
                        source['supplement_url'], source['supplement_sha256'])
        shutil.copy2(archive, root / 'upstream.tar.gz')
        with tarfile.open(archive) as tar:
            for item in tar:
                # Upstream test assets contain absolute author-local symlinks.
                # They are recorded in source.lock.json and never used here.
                if item.issym() or item.islnk():
                    continue
                tar.extract(item, root, filter='data')
        genome = root / 'references'
        genome.mkdir()
        (genome / 'hg38').symlink_to(GENOME, target_is_directory=True)
        (genome / 'default').symlink_to('hg38', target_is_directory=True)
        if sha(GENOME / 'CpG.bed.gz') != CPG_SHA:
            raise ValueError('hg38 CpG namespace changed')
        versions = {'nextflow': run(['nextflow', '-version'], capture_output=True, text=True).stdout,
                    'apptainer': run(['apptainer', '--version'], capture_output=True, text=True).stdout}
        images = read_json(root / 'config/images.lock.json')
        resolved = {}
        for name, specification in images.items():
            target = root / 'images' / (name + '.sif')
            run(['apptainer', 'pull', target, specification['image']], root / 'commands.jsonl')
            resolved[name] = dict(specification, sif=str(target), sha256=sha(target))
        write_json(root / 'images.resolved.json', resolved)
        # Plugin goes into this frozen runtime, never the mutable user plugin cache.
        env = dict(os.environ, NXF_PLUGINS_DIR=str(root / 'plugins'))
        run(['nextflow', 'plugin', 'install', 'nf-schema@2.2.0'], root / 'commands.jsonl', env=env)
        plugin_files = list((root / 'plugins').rglob('*.jar'))
        if not plugin_files:
            raise ValueError('nf-schema plugin installation produced no JAR')
        versions['plugin_files'] = {str(p.relative_to(root)): sha(p) for p in plugin_files}
        bind = ['--bind', f'{root}:{root}', '--bind', f'{PROJECT}:{PROJECT}:ro',
                '--bind', f'{GENOME}:{GENOME}:ro', '--bind', f'{root}/tmp:/tmp']
        env_args = ['--env', f'TMPDIR=/tmp,XDG_CACHE_HOME=/tmp/cache']
        for name, command in {
            'epidish': ['Rscript', '-e', 'write.table(installed.packages()[,c("Package","Version")],stdout(),sep="\\t",row.names=FALSE,quote=FALSE)'],
            'uxm': ['bash', '-c', 'command -v wgbstools; command -v uxm; wgbstools --version; samtools --version | head -n 2'],
            'preprocessing': ['bedtools', '--version'],
            'python': ['python', '--version']
        }.items():
            result = run(['apptainer', 'exec', '--cleanenv', '--no-home', *bind, *env_args,
                          root / 'images' / (name + '.sif'), *command], capture_output=True, text=True)
            (root / (name + '.versions.txt')).write_text(result.stdout + result.stderr)
        # Save exact wrapper source; image labels are not installed-version proof.
        result = run(['apptainer', 'exec', root / 'images/epidish.sif', 'cat', '/source/EpiDISH.R'],
                     capture_output=True, text=True)
        (root / 'EpiDISH.R').write_text(result.stdout)
        write_json(root / 'versions.json', versions)
        # Lock every executable/config resource required by offline execution.
        files = [p for p in root.rglob('*') if p.is_file() and
                 not p.is_relative_to(root / 'tmp') and not p.is_relative_to(genome) and
                 p.name not in ('execution.json', 'runtime.lock.json')]
        write_json(root / 'runtime.lock.json', {str(p.relative_to(root)): sha(p) for p in files})
        finish(root)
        (root / 'STAGED').write_text(utc() + '\n')
    except Exception as exc:
        finish(root, 'failed', error=str(exc))
        raise


if __name__ == '__main__':
    main()
