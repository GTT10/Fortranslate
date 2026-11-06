#!/usr/bin/env python3
"""Copy canonical LiDryer and attach exact pinned TGReact Simple transport fits.

The destination must be new. Chemistry and thermochemistry are preserved;
this preparation does not establish application or physical qualification.
"""
import argparse
import hashlib
import json
from pathlib import Path

from check_mechanism_provenance import check, read_bundle
from pelephysics_transport_import import import_simple_transport, regenerate_simple_transport


def prepare(reference_source: Path, output: Path) -> Path:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / 'references/pelec_tgreact_eb_amr.json').read_text())
    baseline = json.loads((root / 'references/pelec_baseline.json').read_text())
    if manifest['pelec_commit'] != baseline['pelec']['commit']:
        raise ValueError('TGReact and baseline PeleC commits differ')
    canonical = root / 'mechanisms/lidryer.json'
    if not check(canonical):
        raise ValueError('Canonical LiDryer provenance failed')
    bundle = read_bundle(canonical)
    source = root / 'mechanisms' / bundle['source']['file']
    relative = 'Submodules/PelePhysics/Mechanisms/LiDryer/mechanism.yaml'
    native_source = reference_source / relative
    if hashlib.sha256(native_source.read_bytes()).hexdigest() != manifest['source_sha256'][relative]:
        raise ValueError('Pinned LiDryer YAML SHA mismatch')
    if source.read_bytes() != native_source.read_bytes():
        raise ValueError('Canonical and pinned LiDryer YAML differ')
    header = reference_source / manifest['simple_transport_header']['path']
    bundle['simple_transport'] = import_simple_transport(header, bundle['species'], manifest)
    output.mkdir(parents=True, exist_ok=False)
    (output / bundle['source']['file']).write_bytes(source.read_bytes())
    (output / bundle['simple_transport']['source']['file']).write_bytes(header.read_bytes())
    if regenerate_simple_transport(bundle, output) != bundle['simple_transport']:
        raise ValueError('Simple transport round trip differs')
    target = output / 'lidryer.json'
    target.write_text(json.dumps(bundle, indent=2, allow_nan=False) + '\n')
    if not check(target):
        raise ValueError('Prepared bundle provenance failed')
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(prepare(args.reference_source.resolve(), args.output.resolve()))


if __name__ == '__main__':
    main()
