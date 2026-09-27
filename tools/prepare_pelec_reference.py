"""Export pinned PeleC and dependency trees from local Git object databases.

Existing checkouts are read only. The destination must not exist. No network
fetch, checkout, source patching, or build is performed.
"""
import argparse
import io
import json
from pathlib import Path
import subprocess
import tarfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--manifest', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'references/pelec_baseline.json')
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    entries = [{'path': '.', 'commit': manifest['pelec']['commit']}] + manifest['submodules']
    destination = args.destination.resolve()
    if destination.exists():
        parser.error(f'destination already exists: {destination}')
    for entry in entries:
        relative = Path(entry['path'])
        if relative.is_absolute() or '..' in relative.parts:
            parser.error(f'invalid submodule path: {relative}')
        subprocess.run(['git', '-C', str(args.reference / relative), 'cat-file', '-e',
                        entry['commit'] + '^{commit}'], check=True)
    destination.mkdir(parents=True)
    for entry in entries:
        target = destination / entry['path']
        target.mkdir(parents=True, exist_ok=True)
        data = subprocess.check_output(['git', '-C', str(args.reference / entry['path']),
                                        'archive', entry['commit']])
        with tarfile.open(fileobj=io.BytesIO(data)) as stream:
            stream.extractall(target, filter='data')
        print(entry['path'], entry['commit'])
    (destination / 'pelef-reference-provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(destination)


if __name__ == '__main__':
    main()
