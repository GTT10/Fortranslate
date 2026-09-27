"""Run the pinned full HIT application and compare one hydro step in SI units."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('pelec', 'fortran', 'exporter', 'velocity', 'work'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--transport', action='store_true', help='Enable matched constant viscosity and conductivity')
    parser.add_argument('--les', choices=('smagorinsky', 'wale', 'vreman'))
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1]
    for key in ('pelec', 'fortran', 'exporter', 'velocity'):
        setattr(args, key, getattr(args, key).resolve(strict=True))
    velocity_hash = hashlib.sha256(args.velocity.read_bytes()).hexdigest()
    assert velocity_hash == '983ebe8d4db570e56a634346871993e346ba6b3bb7eebf02ca8f39a6e5194776'
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='evolution_', dir=args.work.resolve()))
    print(work, flush=True)
    shutil.copy(args.velocity, work / 'hit_ic_4_32.dat')
    for name in ('pelec_reference.inp', 'hydro_reference.nml'):
        shutil.copy(source / 'cases/reactive_hit_3d' / name, work / name)
    if args.transport:
        reference_input = work / 'pelec_reference.inp'
        reference_input.write_text(reference_input.read_text().replace(
            'pelec.diffuse_vel = 0', 'pelec.diffuse_vel = 1').replace(
            'pelec.diffuse_temp = 0', 'pelec.diffuse_temp = 1'))
        fortran_input = work / 'hydro_reference.nml'
        fortran_input.write_text(fortran_input.read_text().replace(
            'transport_enabled=.false.', 'transport_enabled=.true.') +
            "\n&reactive_constant_transport enabled=.true.,mode='hit',lambda0=.005,reynolds_lambda0=100,prandtl=.71 /\n")
    if args.les:
        model = dict(smagorinsky=0, wale=2, vreman=3)[args.les]
        with (work / 'pelec_reference.inp').open('a') as stream:
            stream.write(f'\npelec.do_les = 1\npelec.les_model = {model}\n'
                         'pelec.Cs = .5\npelec.Cw = .5\npelec.CI = .09\npelec.PrT = .7\n')
        with (work / 'hydro_reference.nml').open('a') as stream:
            stream.write(f"\n&reactive_les enabled=.true.,model='{args.les}',cs=.5,ci=.09,prt=.7 /\n")
    env = dict(os.environ)
    for key in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY'):
        env.pop(key, None)

    def run(command, log):
        with (work / log).open('w') as stream:
            result = subprocess.run([str(v) for v in command], cwd=work, env=env,
                                    stdout=stream, stderr=subprocess.STDOUT, timeout=120)
        assert result.returncode == 0, (work / log).read_text()

    run([args.pelec, 'pelec_reference.inp'], 'pelec.log')
    run([args.fortran, 'hydro_reference.nml'], 'fortran.log')
    run([args.exporter, 'plt00001', 'reference.csv'], 'export.log')
    run([args.exporter, 'plt00000', 'reference_initial.csv'], 'export_initial.log')

    def rows(name):
        with (work / name).open() as stream:
            return list(csv.DictReader(stream))

    reference, actual = rows('reference.csv'), rows('fortran.csv')
    initial = rows('reference_initial.csv')
    assert len(reference) == len(actual) == len(initial) == 64
    gas_r = 8314.46261815324 / 28.97
    rho = 101300 / (gas_r * 300)
    momentum = rho * .1 * math.sqrt(1.4 * gas_r * 300) / math.sqrt(3)
    fields = dict(rho=('density', 1000, rho), rhou=('xmom', 10, momentum),
                  rhov=('ymom', 10, momentum), rhow=('zmom', 10, momentum),
                  rhoE=('rho_E', .1, 101300/.4), temperature=('Temp', 1, 300),
                  pressure=('pressure', .1, 101300))
    for a, b in zip(actual, reference):
        assert float(a['time']) == float(b['time']) == 1e-8
        for key in ('x', 'y', 'z'):
            assert abs(float(a[key]) - .01*float(b[key])) < 1e-16
    errors, updates = {}, {}
    for key, (name, factor, scale) in fields.items():
        errors[key] = max(abs(float(a[key])-factor*float(b[name]))/scale
                          for a, b in zip(actual, reference))
        updates[key] = max(factor*abs(float(b[name])-float(c[name]))/scale
                           for b, c in zip(reference, initial))
    report = dict(transport_enabled=args.transport, les=args.les, normalized_errors=errors, normalized_reference_updates=updates,
                  absolute_normalized_gate=1e-10, velocity_sha256=velocity_hash,
                  binaries={key: hashlib.sha256(getattr(args, key).read_bytes()).hexdigest()
                            for key in ('pelec', 'fortran', 'exporter')})
    (work / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    assert max(updates.values()) > 1e-7, 'Reference did not materially evolve'
    assert all(math.isfinite(v) and v < 1e-10 for v in errors.values()), errors
    print('PASS: full HIT one-step comparison')


if __name__ == '__main__':
    main()
