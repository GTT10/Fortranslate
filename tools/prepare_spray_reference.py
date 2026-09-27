#!/usr/bin/env python3
"""Fetch or verify a pinned, optional fuel mechanism and generate a bundle.

This never writes into mechanisms/ or changes a configured source. The
references are FFCM1_Red methanol and Lu n-dodecane gas mechanisms.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import urllib.error
import urllib.request

from generate_elementary_mechanism import validate
from ingest_cantera_mechanism import build_bundle

ROOT = Path(__file__).resolve().parents[1]


def checked_bytes(data: bytes, entry: dict) -> bytes:
    actual = hashlib.sha256(data).hexdigest()
    if actual != entry['sha256']:
        raise ValueError(f"{entry['path']}: SHA-256 expected {entry['sha256']}, got {actual}")
    if 'git_blob' in entry:
        blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if blob != entry['git_blob']:
            raise ValueError('Source Git blob does not match pinned upstream data')
    return data


def preserve(path: Path, data: bytes) -> None:
    if path.exists():
        if not path.is_file() or path.read_bytes() != data:
            raise ValueError(f'Refusing to replace different existing data: {path}')
    else:
        with path.open('xb') as stream:
            stream.write(data)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fuel', choices=['methanol', 'dodecane', 'heptane-3sp'], default='methanol')
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--yaml', type=Path, help='Offline, exact pinned YAML bytes')
    parser.add_argument('--license', type=Path, help='Offline, exact upstream license bytes')
    args = parser.parse_args()
    if bool(args.yaml) != bool(args.license):
        parser.error('--yaml and --license must be supplied together')
    reference = {'methanol': 'spray_methanol.json', 'dodecane': 'dodecane_lu.json',
                 'heptane-3sp': 'heptane_3sp.json'}[args.fuel]
    manifest = json.loads((ROOT/'references'/reference).read_text())
    prefix = manifest['symbol_prefix']
    data = {}
    for key, entry in manifest['files'].items():
        local = getattr(args, key)
        if local:
            content = local.read_bytes()
        else:
            url = f"https://raw.githubusercontent.com/Pele-Suite/PelePhysics/{manifest['commit']}/{entry['path']}"
            request = urllib.request.Request(url, headers={'User-Agent': 'PeleF-pinned-reference'})
            with urllib.request.urlopen(request, timeout=60) as response:
                content = response.read(2_000_001)
            if len(content)>2_000_000:
                raise ValueError('Unexpectedly large upstream input')
        data[key] = checked_bytes(content, entry)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    for key, entry in manifest['files'].items():
        preserve(out/entry['basename'], data[key])
    bundle = build_bundle(out/manifest['files']['yaml']['basename'], phase=manifest['phase'],
                          module_name=prefix+'_mechanism_mod', symbol_prefix=prefix,
                          chemistry_integrator='implicit', nonreacting=manifest.get('nonreacting', False),
                          description=('Pinned FFCM1_Red: optional methanol spray verification, not diesel validation'
                                       if args.fuel == 'methanol' else
                                       'Pinned reaction-free NC7H16/O2/N2 spray reference' if args.fuel == 'heptane-3sp' else
                                       'Pinned Lu n-dodecane: gas chemistry integration, not experimental spray validation'))
    validate(bundle)
    if len(bundle['species']) != manifest['species_count'] or len(bundle['reactions']) != manifest['reaction_count']:
        raise ValueError('Unexpected normalized mechanism dimensions')
    preserve(out/(prefix+'.json'), (json.dumps(bundle, indent=2, sort_keys=True, allow_nan=False)+'\n').encode())
    preserve(out/'upstream.json', (json.dumps(manifest, indent=2)+'\n').encode())
    print(f"PASS: {out/(prefix+'.json')} ({manifest['species_count']} species, "
          f"{manifest['reaction_count']} reactions; source and license verified)")
    return 0


if __name__=='__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, urllib.error.URLError) as error:
        print(f'prepare_spray_reference: {error}', file=sys.stderr)
        raise SystemExit(1)
