#!/usr/bin/env python3
"""Check coupled wall/outflow budgets and physical-boundary checkpoint replay."""
import argparse
import csv
import re
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--inflow', action='store_true')
    parser.add_argument('--sgs', type=float, default=0.0)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='spray-boundaries-', dir=args.work_dir.resolve()))
    case = 'inflow.nml' if args.inflow else 'wall_reflection.nml'
    template = (Path(__file__).resolve().parents[1] / 'cases/spray_3d' / case).read_text()
    template = template.replace('smagorinsky_constant=0,', f'smagorinsky_constant={args.sgs},')
    if args.sgs:
        template = template.replace('&spray', '&spray\n shear_velocity=10,')
    for kind in (('inflow',) if args.inflow else ('wall', 'outflow')):
        def run(name, extra=''):
            tag = kind + '_' + name
            config = re.sub(r"output_prefix='[^']*'", f"output_prefix='{tag}'", template)
            if kind == 'outflow':
                config = config.replace("'slip_wall'", "'outflow'")
            config = config.replace('&spray', '&spray\n' + extra)
            path = work / (tag + '.nml')
            path.write_text(config)
            result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                    text=True, capture_output=True, timeout=180)
            (work / (tag + '.log')).write_text(result.stdout + result.stderr)
            if result.returncode:
                raise RuntimeError(result.stdout + result.stderr)
            return tag
        ref = run('reference')
        run('stop', f"stop_after_steps=2, checkpoint_file='{kind}.chk',")
        resumed = run('resumed', f"restart_file='{kind}.chk',")
        for suffix in ('.gas.csv', '.parcels.csv'):
            assert (work / (ref + suffix)).read_bytes() == (work / (resumed + suffix)).read_bytes(), suffix
        def last(tag):
            with (work / (tag + '.history.csv')).open() as stream:
                return list(csv.DictReader(stream))[-1]
        row = last(ref)
        assert row == last(resumed), 'resumed history ledger'
        mass = float(row['system_mass']) + abs(float(row['boundary_mass']))
        energy = abs(float(row['total_energy'])) + abs(float(row['boundary_energy']))
        assert abs(float(row['mass_error'])) < 1e-9 * mass
        assert abs(float(row['energy_error'])) < 1e-9 * energy
        for axis in ('px', 'py', 'pz'):
            assert abs(float(row[axis + '_error'])) < 1e-9 * mass * 100
        if kind == 'inflow':
            assert float(row['boundary_mass']) > 0
            with (work / (ref + '.gas.csv')).open() as stream:
                assert sum(float(cell['rhoY_H2']) for cell in csv.DictReader(stream)) > 0
        elif kind == 'outflow':
            assert int(row['live_parcels']) == 0
            assert float(row['boundary_mass']) < -5e-11
        else:
            assert int(row['live_parcels']) == 2
            assert float(row['boundary_pz']) < -1e-9
        if kind in ('wall', 'inflow'):
            bad_input = (work / (resumed + '.nml')).read_text()
            bad_input = (bad_input.replace("'slip_wall'", "'outflow'") if kind == "wall"
                         else bad_input.replace("inflow_temperature(1)=1100", "inflow_temperature(1)=1101"))
            bad_input = bad_input.replace(resumed, 'changed_faces')
            bad_path = work / 'changed_faces.nml'
            bad_path.write_text(bad_input)
            result = subprocess.run([str(args.executable.resolve()), str(bad_path)], cwd=work,
                                    text=True, capture_output=True, timeout=30)
            (work / 'changed_faces.log').write_text(result.stdout + result.stderr)
            assert result.returncode != 0 and 'incompatible spray checkpoint' in result.stdout + result.stderr
            assert not (work / 'changed_faces.history.csv').exists()
        print(kind, 'PASS: exact restart and boundary budget')
    print('artifacts:', work)


if __name__ == '__main__':
    main()
