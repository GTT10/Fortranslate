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
    parser.add_argument('--ppm', action='store_true')
    parser.add_argument('--difmag', type=float)
    parser.add_argument('--mach-inlet', action='store_true')
    parser.add_argument('--thermal-wall', action='store_true')
    parser.add_argument('--acceleration', type=float, default=0.0)
    parser.add_argument('--sgs', type=float, default=0.0)
    parser.add_argument('--pressure-outlet', action='store_true')
    parser.add_argument('--les-isotropic', action='store_true')
    parser.add_argument('--les-harmonic', action='store_true')
    parser.add_argument('--dynamic-physical', action='store_true')
    parser.add_argument('--les-model', choices=('wale', 'vreman', 'dynamic'))
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='spray-boundaries-', dir=args.work_dir.resolve()))
    case = 'thermal_moving_wall.nml' if args.thermal_wall else ('inflow.nml' if args.inflow else 'wall_reflection.nml')
    template = (Path(__file__).resolve().parents[1] / 'cases/spray_3d' / case).read_text()
    if args.ppm:
        template = template.replace('&spray', "&spray\n hydro_reconstruction='characteristic_ppm',")
    if args.difmag is not None:
        if args.difmag < 0:
            parser.error('--difmag requires a nonnegative coefficient')
        template = template.replace('&spray', f"&spray\n difmag={args.difmag}, riemann_solver='pelec',")
    template = template.replace('smagorinsky_constant=0,', f'smagorinsky_constant={args.sgs},')
    if args.les_isotropic:
        if args.les_model not in ('wale', 'vreman'):
            parser.error('--les-isotropic requires --les-model wale or vreman')
        template = template.replace('&spray', '&spray\n les_ci=0.1,')
    if args.les_harmonic:
        if args.les_model != 'dynamic':
            parser.error('--les-harmonic requires --les-model dynamic')
        template = template.replace('&spray', '&spray\n les_harmonic_mean=.true.,')
    if args.les_model:
        template = template.replace('&spray', f"&spray\n les_model='{args.les_model}',")
    if args.dynamic_physical and args.les_model != 'dynamic':
        parser.error('--dynamic-physical requires --les-model dynamic')
    if args.les_model == 'dynamic':
        if not args.dynamic_physical:
            template = template.replace("'slip_wall'", "'periodic'")
        template = template.replace('&spray', "&spray\n les_test_filter='gaussian', les_test_ratio=4, les_grid_ratio=2,")
    if args.mach_inlet:
        if not args.inflow or args.thermal_wall:
            parser.error('--mach-inlet requires --inflow without --thermal-wall')
        if args.les_model == 'dynamic' and not args.dynamic_physical:
            parser.error('--mach-inlet with dynamic requires --dynamic-physical')
        template = template.replace('inflow_velocity(:,1)=200,0,0', 'inflow_velocity(:,1)=0,0,0')
        template = template.replace('inflow_pressure(1)=101325', 'inflow_pressure(1)=0')
        template = template.replace('&spray', '&spray\n inflow_mach(1)=0.2,')
    if args.acceleration:
        template = template.replace('&spray', f'&spray\n acceleration(3)={args.acceleration},')
    if args.sgs:
        template = template.replace('&spray', '&spray\n shear_velocity=10,')
    kinds = ('thermal_wall',) if args.thermal_wall else (('inflow',) if args.inflow else ('wall', 'outflow'))
    if args.les_model == 'dynamic' and not args.dynamic_physical:
        kinds = ('periodic',)
    if args.pressure_outlet:
        if (args.les_model == 'dynamic' and not args.dynamic_physical) or args.inflow or args.thermal_wall:
            parser.error('--pressure-outlet requires a physical outflow fixture')
        kinds = ('outflow',)
        template = template.replace('&spray', '&spray\n outlet_pressure(6)=100000,')
        template = template.replace('transport_enabled=.false.', 'transport_enabled=.true.')
    for kind in kinds:
        def run(name, extra=''):
            tag = kind + '_' + name
            config = re.sub(r"output_prefix='[^']*'", f"output_prefix='{tag}'", template)
            if kind == 'outflow':
                config = config.replace("'slip_wall'", "'outflow'")
            config = config.rsplit('/', 1)[0] + extra + '\n/\n'
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
        if args.acceleration:
            assert float(row['source_pz']) < 0
        boundary_mass = 0 if kind == 'periodic' else float(row['boundary_mass'])
        boundary_energy = 0 if kind == 'periodic' else float(row['boundary_energy'])
        mass = float(row['system_mass']) + abs(boundary_mass)
        energy = abs(float(row['total_energy'])) + abs(boundary_energy)
        assert abs(float(row['mass_error'])) < 1e-9 * mass
        assert abs(float(row['energy_error'])) < 1e-9 * energy
        for axis in ('px', 'py', 'pz'):
            assert abs(float(row[axis + '_error'])) < 1e-9 * mass * 100
        if kind == 'inflow':
            assert float(row['boundary_mass']) > 0
            with (work / (ref + '.gas.csv')).open() as stream:
                h2 = sum(float(cell['rhoY_H2']) for cell in csv.DictReader(stream))
                if args.mach_inlet:
                    assert h2 == 0, 'Mach inlet uses interior composition, not the H2 reservoir'
                else:
                    assert h2 > 0
        elif kind == 'thermal_wall':
            assert float(row['boundary_energy']) > 0
            assert int(row['live_parcels']) == 2
        elif kind == 'outflow':
            assert int(row['live_parcels']) == 0
            assert float(row['boundary_mass']) < -5e-11
        elif kind == 'periodic':
            assert int(row['live_parcels']) == 2
            assert not any(name.startswith('boundary_') for name in row)
        else:
            assert int(row['live_parcels']) == 2
            assert float(row['boundary_pz']) < -1e-9
        if kind in ('wall', 'inflow', 'thermal_wall'):
            bad_input = (work / (resumed + '.nml')).read_text()
            if kind == 'thermal_wall':
                bad_input = bad_input.replace('wall_velocity(1,6)=10', 'wall_velocity(1,6)=11')
            bad_input = (bad_input.replace("'slip_wall'", "'outflow'") if kind == "wall"
                         else bad_input.replace("inflow_temperature(1)=1100", "inflow_temperature(1)=1101"))
            if args.acceleration:
                bad_input = (work / (resumed + '.nml')).read_text().replace(
                    f'acceleration(3)={args.acceleration}', f'acceleration(3)={args.acceleration + 1}')
            bad_input = bad_input.replace(resumed, 'changed_faces')
            bad_path = work / 'changed_faces.nml'
            bad_path.write_text(bad_input)
            result = subprocess.run([str(args.executable.resolve()), str(bad_path)], cwd=work,
                                    text=True, capture_output=True, timeout=30)
            (work / 'changed_faces.log').write_text(result.stdout + result.stderr)
            assert result.returncode != 0 and 'incompatible spray checkpoint' in result.stdout + result.stderr
            assert not (work / 'changed_faces.history.csv').exists()
        if args.ppm:
            changed = (work / (resumed + '.nml')).read_text().replace(
                "hydro_reconstruction='characteristic_ppm'", "hydro_reconstruction='characteristic_plm'")
            changed = changed.replace(resumed, 'changed_reconstruction')
            path = work / 'changed_reconstruction.nml'
            path.write_text(changed)
            result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                    text=True, capture_output=True, timeout=30)
            (work / 'changed_reconstruction.log').write_text(result.stdout + result.stderr)
            diagnostic = 'incompatible spray checkpoint'
            assert result.returncode != 0 and diagnostic in result.stdout + result.stderr
            assert not (work / 'changed_reconstruction.history.csv').exists()
        if args.difmag is not None:
            original = (work / (resumed + '.nml')).read_text()
            for label, replacement, diagnostic in [
                ('changed_difmag', f'difmag={args.difmag + 0.1},', 'incompatible spray checkpoint'),
                ('invalid_difmag', 'difmag=-2,', 'Invalid artificial-viscosity coefficient'),
                ('nonfinite_difmag', 'difmag=NaN,', 'Invalid artificial-viscosity coefficient'),
            ]:
                config = original.replace(f'difmag={args.difmag},', replacement).replace(resumed, label)
                path = work / (label + '.nml')
                path.write_text(config)
                result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                        text=True, capture_output=True, timeout=30)
                (work / (label + '.log')).write_text(result.stdout + result.stderr)
                assert result.returncode != 0 and diagnostic in result.stdout + result.stderr
                assert not (work / (label + '.history.csv')).exists()
            # Exercise the disabled correction through the complete coupled path.
            zero = run('zero_difmag', 'difmag=0,')
            assert (work / (zero + '.gas.csv')).read_bytes() != (work / (ref + '.gas.csv')).read_bytes()
            # Equivalent automatic and explicit default settings must share restart identity.
            run('auto_stop', f"difmag=-1, stop_after_steps=2, checkpoint_file='{kind}_auto.chk',")
            auto = run('auto_resume', f"difmag=-1, restart_file='{kind}_auto.chk',")
            explicit = run('explicit_resume', f"difmag={0.1 if args.ppm else 0}, restart_file='{kind}_auto.chk',")
            for suffix in ('.gas.csv', '.parcels.csv', '.history.csv'):
                assert (work / (auto + suffix)).read_bytes() == (work / (explicit + suffix)).read_bytes()
        if args.mach_inlet:
            original = (work / (resumed + '.nml')).read_text()
            for label, changed, diagnostic in [
                ('changed_mach', original.replace('inflow_mach(1)=0.2', 'inflow_mach(1)=0.3'),
                 'incompatible spray checkpoint'),
                ('changed_mach_pressure', original.replace('inflow_pressure(1)=0', 'inflow_pressure(1)=101325'),
                 'incompatible spray checkpoint'),
                ('invalid_mach', original.replace('inflow_mach(1)=0.2', 'inflow_mach(1)=101'),
                 'Active inlet Mach number'),
            ]:
                changed = changed.replace(resumed, label)
                path = work / (label + '.nml')
                path.write_text(changed)
                result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                        text=True, capture_output=True, timeout=30)
                (work / (label + '.log')).write_text(result.stdout + result.stderr)
                assert result.returncode != 0 and diagnostic in result.stdout + result.stderr
                assert not (work / (label + '.history.csv')).exists()
        if args.pressure_outlet:
            original = (work / (resumed + '.nml')).read_text()
            for label, replacement, diagnostic in [
                ('changed_pressure', 'outlet_pressure(6)=100100,', 'incompatible spray checkpoint'),
                ('removed_pressure', '', 'incompatible spray checkpoint'),
                ('invalid_pressure', 'outlet_pressure(6)=-1,', 'Invalid outlet pressure'),
            ]:
                config = original.replace('outlet_pressure(6)=100000,', replacement).replace(resumed, label)
                path = work / (label + '.nml')
                path.write_text(config)
                result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                        text=True, capture_output=True, timeout=30)
                (work / (label + '.log')).write_text(result.stdout + result.stderr)
                assert result.returncode != 0 and diagnostic in result.stdout + result.stderr
                assert not (work / (label + '.history.csv')).exists()
        if args.les_model:
            original = (work / (resumed + '.nml')).read_text()
            for label, replacement, diagnostic in [
                ('changed_model', "les_model='smagorinsky',", 'incompatible spray checkpoint'),
                ('removed_model', '', 'incompatible spray checkpoint'),
                ('invalid_model', "les_model='invalid',", 'Invalid LES model'),
            ]:
                config = original.replace(f"les_model='{args.les_model}',", replacement).replace(resumed, label)
                if args.les_model == 'dynamic' and label != 'invalid_model':
                    config = config.replace("les_test_filter='gaussian', les_test_ratio=4, les_grid_ratio=2,", '')
                    config = config.replace('les_harmonic_mean=.true.,', '')
                path = work / (label + '.nml')
                path.write_text(config)
                result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                        text=True, capture_output=True, timeout=30)
                (work / (label + '.log')).write_text(result.stdout + result.stderr)
                assert result.returncode != 0 and diagnostic in result.stdout + result.stderr
                assert not (work / (label + '.history.csv')).exists()
        if args.les_isotropic:
            for label, replacement, diagnostic in [
                ('changed_ci', 'les_ci=0.2,', 'incompatible spray checkpoint'),
                ('removed_ci', '', 'incompatible spray checkpoint'),
                ('invalid_ci', 'les_ci=-0.1,', 'Invalid LES CI'),
            ]:
                config = original.replace('les_ci=0.1,', replacement).replace(resumed, label)
                path = work / (label + '.nml')
                path.write_text(config)
                result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                        text=True, capture_output=True, timeout=30)
                (work / (label + '.log')).write_text(result.stdout + result.stderr)
                assert result.returncode != 0 and diagnostic in result.stdout + result.stderr
                assert not (work / (label + '.history.csv')).exists()
        if args.les_model == 'dynamic':
            guards = [
                ('filter', "les_test_filter='gaussian'", "les_test_filter='box'", 'incompatible spray checkpoint'),
                ('ratio', 'les_test_ratio=4', 'les_test_ratio=2', 'incompatible spray checkpoint'),
                ('grid', 'les_grid_ratio=2', 'les_grid_ratio=1', 'incompatible spray checkpoint'),
                ('removed_filter', "les_test_filter='gaussian',", '', 'incompatible spray checkpoint'),
                ('removed_ratio', 'les_test_ratio=4,', '', 'incompatible spray checkpoint'),
                ('removed_grid', 'les_grid_ratio=2,', '', 'incompatible spray checkpoint'),
                ('bad_filter', 'les_test_ratio=4', 'les_test_ratio=3', 'Invalid dynamic LES filter'),
                ('boundary', "'periodic'", "'outflow'", 'incompatible spray checkpoint'),
            ]
            if args.dynamic_physical:
                guards = [guard for guard in guards if guard[0] != 'boundary']
            if args.les_harmonic:
                guards += [
                    ('harmonic', 'les_harmonic_mean=.true.', 'les_harmonic_mean=.false.', 'incompatible spray checkpoint'),
                    ('removed_harmonic', 'les_harmonic_mean=.true.,', '', 'incompatible spray checkpoint'),
                ]
            for label, old, new, diagnostic in guards:
                config = original.replace(old, new).replace(resumed, label)
                path = work / (label + '.nml')
                path.write_text(config)
                result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                        text=True, capture_output=True, timeout=30)
                (work / (label + '.log')).write_text(result.stdout + result.stderr)
                assert result.returncode != 0 and diagnostic in result.stdout + result.stderr
                assert not (work / (label + '.history.csv')).exists()
        print(kind, 'PASS: exact restart and boundary budget')
    print('artifacts:', work)


if __name__ == '__main__':
    main()
