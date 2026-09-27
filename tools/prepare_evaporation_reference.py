"""Prepare an in-range 3D adaptation of the pinned Spray-Evaporation reference.

This prepares reference data only; it does not claim trajectory parity with
PeleF's serial spray integrator or reproduce upstream's original 273 K case.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    bundle = json.loads(args.bundle.read_text())
    manifest = json.loads((root / 'references/heptane_3sp.json').read_text())
    assert bundle.get('nonreacting') is True and bundle['reactions'] == []
    assert bundle['species'] == ['NC7H16', 'O2', 'N2']
    assert bundle['source']['sha256'] == manifest['files']['yaml']['sha256']
    source = args.bundle.parent / bundle['source']['file']
    assert hashlib.sha256(source.read_bytes()).hexdigest() == bundle['source']['sha256']
    mw = {s['name']: s['molecular_weight'] for s in bundle['thermo']}
    denominator = .21*mw['O2'] + .79*mw['N2']
    oxygen = .21*mw['O2']/denominator
    multiplicity = 1e-10 / (4 * 681.4 * math.pi/6 * (20e-6)**3)
    args.output.mkdir(parents=True, exist_ok=False)
    text = (root / 'cases/spray_reference/evaporation.inp.in').read_text()
    text = text.replace('@O2@', format(oxygen, '.17g')).replace('@N2@', format(1-oxygen, '.17g'))
    (args.output / 'reference.inp').write_text(text)
    # 3 positions + 11 SprayComps real values: velocity, T, diameter,
    # composition, current/initial multiplicity, breakup state and film height.
    row = [.15,.15,.15,0,0,200,310,.002,1,multiplicity,multiplicity,0,0,0]
    (args.output / 'particles.dat').write_text('4\n' + (' '.join(format(v, '.17g') for v in row)+'\n')*4)
    (args.output / 'preparation.json').write_text(json.dumps(dict(
        reference=manifest, gas_temperature=400, liquid_temperature=310,
        total_liquid_mass_kg=1e-10, multiplicity=multiplicity,
        gas_mass_fractions=dict(O2=oxygen,N2=1-oxygen),
        note='3D periodic in-range adaptation; time-integrator parity not established'), indent=2)+'\n')
    print(args.output.resolve())


if __name__ == '__main__':
    main()
