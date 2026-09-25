#!/usr/bin/env python3
"""Focused input-to-output check for the fixed MPI spray executable."""
import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith('#')))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--mpiexec', type=Path, required=True)
    parser.add_argument('--ranks', type=int, required=True)
    parser.add_argument('--species-count', type=int, default=10)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--pressure-outlet', action='store_true')
    parser.add_argument('--rate-files', action='store_true')
    parser.add_argument('--ppm', action='store_true')
    parser.add_argument('--difmag', type=float)
    parser.add_argument('--les-isotropic', action='store_true')
    parser.add_argument('--les-harmonic', action='store_true')
    parser.add_argument('--dynamic-walls', action='store_true')
    parser.add_argument('--mach-inlet', action='store_true')
    parser.add_argument('--les-model', choices=('wale', 'vreman', 'dynamic'))
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='run-', dir=args.work_dir.resolve()))
    source = args.input.read_text()
    if args.ppm:
        source = source.replace('&mpi_spray', "&mpi_spray\n hydro_reconstruction='characteristic_ppm',")
    if args.difmag is not None:
        if not math.isfinite(args.difmag) or args.difmag < 0:
            parser.error('--difmag requires a finite nonnegative coefficient')
        source = source.replace('&mpi_spray', f"&mpi_spray\n difmag={args.difmag}, solver='pelec',")
    input_path = args.input.resolve()
    if args.pressure_outlet:
        if args.les_model not in ('wale', 'vreman'):
            parser.error('--pressure-outlet requires a static LES fixture')
        source = source.replace('&mpi_spray', "&mpi_spray\n boundary_faces=6*'outflow', outlet_pressure(2)=320000,")
    if args.les_isotropic:
        if args.les_model not in ('wale', 'vreman'):
            parser.error('--les-isotropic requires --les-model wale or vreman')
        source = source.replace('&mpi_spray', '&mpi_spray\n les_ci=0.1,')
    if args.les_harmonic:
        if args.les_model != 'dynamic':
            parser.error('--les-harmonic requires --les-model dynamic')
        source = source.replace('&mpi_spray', '&mpi_spray\n les_harmonic_mean=.true.,')
    if args.dynamic_walls:
        if args.les_model != 'dynamic':
            parser.error('--dynamic-walls requires --les-model dynamic')
        source = source.replace('&mpi_spray', "&mpi_spray\n boundary_faces='periodic','periodic',"
                                "'no_slip_wall','isothermal_wall','periodic','periodic',"
                                "wall_temperature=6*600.,wall_velocity(1,3)=2.,")
    if args.mach_inlet:
        if args.les_model != 'dynamic' or args.dynamic_walls or args.pressure_outlet:
            parser.error('--mach-inlet requires dynamic without other boundary fixtures')
        source = source.replace('&mpi_spray', "&mpi_spray\n boundary_faces='inflow','outflow',"
                                "'periodic','periodic','periodic','periodic',"
                                "inflow_mach(1)=0.2,inflow_temperature(1)=900,inflow_pressure(1)=0,")
    if args.les_model:
        source = source.replace('&mpi_spray', f"&mpi_spray\n les_cs=.17, les_model='{args.les_model}',")
        if args.les_model == 'dynamic':
            source = source.replace('&mpi_spray', "&mpi_spray\n les_test_filter='gaussian', les_test_ratio=4, les_grid_ratio=2,")
        input_path = work / 'les.nml'
        input_path.write_text(source)
    if args.ppm or args.difmag is not None:
        input_path = work / 'hydro.nml'
        input_path.write_text(source)
    command = [str(args.mpiexec.resolve()), '-n', str(args.ranks),
               str(args.executable.resolve()), str(input_path)]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=90)
    (work / 'run.log').write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    history = rows(work / 'periodic_mixture.history.csv')
    final_step = int(history[-1]['step'])
    assert final_step >= 2 and len(history) == final_step + 1
    if args.difmag is None:
        assert final_step == 2
    split_time = history[1]['time'] if args.difmag is not None else '1e-8'
    assert math.isclose(float(history[-1]['time']), 2e-8, rel_tol=1e-14)
    for row in history:
        assert all(math.isfinite(float(value)) for value in row.values())
        assert abs(float(row['relative_mass_error'])) < 2e-11
        assert abs(float(row['relative_energy_error'])) < 2e-11
        assert int(row['parcels']) == 4
    if args.mach_inlet:
        assert float(history[-1]['boundary_mass']) > 0
    if args.dynamic_walls:
        assert abs(float(history[-1]['boundary_energy'])) > 1e-14
    gas = [row for path in sorted(work.glob('*.gas.csv')) for row in rows(path)]
    assert len(gas) == 8 * 4 * 4
    cells = {(int(row['i']), int(row['j']), int(row['k'])) for row in gas}
    assert cells == {(i, j, k) for i in range(1, 9) for j in range(1, 5) for k in range(1, 5)}
    for row in gas:
        assert all(math.isfinite(float(value)) for value in row.values())
        rho = float(row['rho'])
        species = [float(value) for name, value in row.items() if name.startswith('rhoY_')]
        assert len(species) == args.species_count
        assert rho > 0 and float(row['temperature']) > 0 and min(species) >= 0
        assert math.isclose(sum(species), rho, rel_tol=2e-12)
    particles = [row for path in sorted(work.glob('*.particles.csv')) for row in rows(path)]
    assert sorted(int(row['id']) for row in particles) == [1, 2, 3, 4]
    for row in particles:
        assert all(math.isfinite(float(value)) for value in row.values())
        assert math.isclose(float(row['mass']), float(row['component_mass_1']) + float(row['component_mass_2']), rel_tol=2e-12)
        assert float(row['temperature']) > 330
    # A second invocation must fail before replacing the saved history.
    original = (work / 'periodic_mixture.history.csv').read_bytes()
    repeated = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=90)
    (work / 'repeat.log').write_text(repeated.stdout + repeated.stderr)
    assert repeated.returncode != 0
    assert (work / 'periodic_mixture.history.csv').read_bytes() == original

    def run_variant(name, rank_count, additions, stop=None, expect_success=True, file_group=''):
        text = source.replace("output_prefix='periodic_mixture'", f"output_prefix='{name}'")
        if stop is not None:
            text = text.replace('final_time=2e-8', f'final_time={stop}')
        text = text.replace('\n/', '\n' + additions + '\n/')
        path = work / f'{name}.nml'
        path.write_text(text + file_group)
        invocation = [str(args.mpiexec.resolve()), '-n', str(rank_count),
                      str(args.executable.resolve()), str(path)]
        output = subprocess.run(invocation, cwd=work, capture_output=True, text=True, timeout=90)
        (work / f'{name}.log').write_text(output.stdout + output.stderr)
        assert (output.returncode == 0) == expect_success, output.stdout + output.stderr

    def snapshots(prefix, suffix, keys):
        records = [row for path in sorted(work.glob(f'{prefix}.rank*.{suffix}.csv')) for row in rows(path)]
        return {tuple(int(row[key]) for key in keys): row for row in records}

    if args.rate_files:
        ports = ('port_count=2,port_events=2,2,port_start=0,7.5e-9,'
                 'port_interval=2*7.5e-9,port_number=2,3,')
        inline = (ports + 'injection_profile_points=2,injection_profile_times=0,1.5e-8,'
                  'injection_profile_rates=.01,.01,port_profile_points=0,3,'
                  'port_profile_times(:,2)=7.5e-9,1.5e-8,2.25e-8,port_profile_rates(:,2)=0,0,.04,')
        (work / 'shared.dat').write_text('# common\n0 .01\n1.5e-8 .01\n')
        (work / 'port2.dat').write_text('7.5e-9 0\n1.5e-8\t0\n2.25e-8 .04 ! end\n')
        group = "\n&mpi_spray_rate_files profile_file='shared.dat',profile_files(2)='port2.dat' /\n"
        run_variant('inline', args.ranks, inline)
        run_variant('files', args.ranks, ports, file_group=group)
        run_variant('split_files', args.ranks, ports + "checkpoint_path='files.chk'",
                    stop='7.5e-9', file_group=group)
        (work / 'renamed.dat').write_text('0 .01\n1.5e-8 .01 # renamed\n')
        run_variant('resume_files', 6 - args.ranks, ports + "restart_path='files.chk'",
                    file_group=group.replace('shared.dat', 'renamed.dat'))
        for name in ('files', 'resume_files'):
            for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
                reference, actual = snapshots('inline', suffix, keys), snapshots(name, suffix, keys)
                assert reference.keys() == actual.keys()
                for key, record in reference.items():
                    for column, value in record.items():
                        assert math.isclose(float(value), float(actual[key][column]), rel_tol=2e-11, abs_tol=1e-25)
        history = rows(work / 'files.history.csv')
        assert math.isclose(float(history[-1]['injected_mass']), 3e-10, rel_tol=1e-14)
        (work / 'shared.dat').write_text('0 .01\n1.5e-8 .02\n')
        run_variant('changed', args.ranks, ports + "restart_path='files.chk'", expect_success=False, file_group=group)
        assert 'Restart rejected' in (work / 'changed.log').read_text()
        (work / 'shared.dat').write_text('0 .01\n1.5e-8 -.01\n')
        run_variant('invalid_file', args.ranks, ports, expect_success=False, file_group=group)
        assert 'Invalid spray rate file' in (work / 'invalid_file.log').read_text()
        for name in ('changed', 'invalid_file'):
            assert not (work / f'{name}.history.csv').exists()
        print(f'PASS MPI rate files, inline parity and changed-rank restart: {work}')
        return

    run_variant('split', args.ranks, " checkpoint_path='split.chk'", stop=split_time)
    assert (work / 'split.chk').is_file()
    reference_gas = snapshots('periodic_mixture', 'gas', ['i', 'j', 'k'])
    reference_particles = snapshots('periodic_mixture', 'particles', ['id'])
    for name, count in [('resumed', args.ranks), ('redistributed', 6 - args.ranks)]:
        run_variant(name, count, " restart_path='split.chk'")
        for reference, actual in [
            (reference_gas, snapshots(name, 'gas', ['i', 'j', 'k'])),
            (reference_particles, snapshots(name, 'particles', ['id'])),
        ]:
            assert actual.keys() == reference.keys()
            for key, expected in reference.items():
                for column, value in expected.items():
                    if name == 'resumed':
                        assert actual[key][column] == value, (name, key, column)
                    else:
                        assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (name, key, column)
        resumed_history = rows(work / f'{name}.history.csv')
        assert int(resumed_history[0]['step']) == 1 and int(resumed_history[-1]['step']) == final_step
        assert int(resumed_history[-1]['parcels']) == 4
        assert abs(float(resumed_history[-1]['relative_energy_error'])) < 2e-11
    if args.mach_inlet:
        original_source = source
        for label, old, new in [
            ('mach', 'inflow_mach(1)=0.2', 'inflow_mach(1)=0.3'),
            ('mach_pressure', 'inflow_pressure(1)=0', 'inflow_pressure(1)=101325'),
            ('mach_temperature', 'inflow_temperature(1)=900', 'inflow_temperature(1)=901'),
            ('mach_removed', 'inflow_mach(1)=0.2,', ''),
        ]:
            source = original_source.replace(old, new)
            run_variant(label, args.ranks, " restart_path='split.chk'", expect_success=False)
            assert 'Restart rejected' in (work / f'{label}.log').read_text()
            assert not (work / f'{label}.history.csv').exists()
        source = original_source
    if args.ppm:
        run_variant('changed_reconstruction', args.ranks,
                    " restart_path='split.chk', hydro_reconstruction='characteristic_plm'", expect_success=False)
        diagnostic = 'Restart rejected'
        assert diagnostic in (work / 'changed_reconstruction.log').read_text()
    run_variant('incompatible', args.ranks, " restart_path='split.chk'\n liquid_cp=4000,2000", expect_success=False)
    assert not (work / 'incompatible.history.csv').exists()
    if args.difmag is not None:
        for label, value, diagnostic in [
            ('changed_difmag', args.difmag + .1, 'Restart rejected'),
            ('invalid_difmag', -2, 'Invalid artificial-viscosity coefficient'),
            ('nonfinite_difmag', 'NaN', 'Invalid artificial-viscosity coefficient'),
        ]:
            run_variant(label, args.ranks, f"restart_path='split.chk', difmag={value}", expect_success=False)
            assert diagnostic in (work / f'{label}.log').read_text()
            assert not (work / f'{label}.history.csv').exists()
        run_variant('auto_split', args.ranks, "difmag=-1, checkpoint_path='auto.chk'", stop='1e-8')
        run_variant('auto_resume', 6 - args.ranks, "difmag=-1, restart_path='auto.chk'")
        run_variant('explicit_resume', 6 - args.ranks, f"difmag={.1 if args.ppm else 0}, restart_path='auto.chk'")
        for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
            assert snapshots('auto_resume', suffix, keys) == snapshots('explicit_resume', suffix, keys)
        assert rows(work / 'auto_resume.history.csv') == rows(work / 'explicit_resume.history.csv')
        print('PASS MPI difmag input and changed-rank restart; artifacts:', work)
        return
    if args.les_model:
        original_source = source
        for label, replacement, diagnostic in [
            ('changed_model', "les_model='smagorinsky',", 'Restart rejected: invalid data or incompatible model/input context'),
            ('removed_model', '', 'Restart rejected: invalid data or incompatible model/input context'),
            ('invalid_model', "les_model='invalid',", 'Invalid LES model'),
        ]:
            source = original_source.replace(f"les_model='{args.les_model}',", replacement)
            if args.les_model == 'dynamic' and label != 'invalid_model':
                source = source.replace("les_test_filter='gaussian', les_test_ratio=4, les_grid_ratio=2,", '')
                source = source.replace('les_harmonic_mean=.true.,', '')
            run_variant(label, args.ranks, " restart_path='split.chk'", expect_success=False)
            assert diagnostic in (work / f'{label}.log').read_text()
            assert not (work / f'{label}.history.csv').exists()
        if args.pressure_outlet:
            for label, replacement, diagnostic in [
                ('changed_pressure', 'outlet_pressure(2)=320100,', 'Restart rejected: invalid data or incompatible model/input context'),
                ('removed_pressure', '', 'Restart rejected: invalid data or incompatible model/input context'),
                ('invalid_pressure', 'outlet_pressure(2)=-1,', 'Invalid outlet pressure'),
            ]:
                source = original_source.replace('outlet_pressure(2)=320000,', replacement)
                run_variant(label, args.ranks, " restart_path='split.chk'", expect_success=False)
                assert diagnostic in (work / f'{label}.log').read_text()
                assert not (work / f'{label}.history.csv').exists()
        if args.les_isotropic:
            for label, replacement, diagnostic in [
                ('changed_ci', 'les_ci=0.2,', 'Restart rejected: invalid data or incompatible model/input context'),
                ('removed_ci', '', 'Restart rejected: invalid data or incompatible model/input context'),
                ('invalid_ci', 'les_ci=-0.1,', 'Invalid LES CI'),
            ]:
                source = original_source.replace('les_ci=0.1,', replacement)
                run_variant(label, args.ranks, " restart_path='split.chk'", expect_success=False)
                assert diagnostic in (work / f'{label}.log').read_text()
                assert not (work / f'{label}.history.csv').exists()
        if args.les_model == 'dynamic':
            guards = [
                ('filter', "les_test_filter='gaussian'", "les_test_filter='box'"),
                ('ratio', 'les_test_ratio=4', 'les_test_ratio=2'),
                ('grid', 'les_grid_ratio=2', 'les_grid_ratio=1'),
                ('removed_filter', "les_test_filter='gaussian',", ''),
                ('removed_ratio', 'les_test_ratio=4,', ''),
                ('removed_grid', 'les_grid_ratio=2,', ''),
            ]
            if args.les_harmonic:
                guards += [
                    ('harmonic', 'les_harmonic_mean=.true.', 'les_harmonic_mean=.false.'),
                    ('removed_harmonic', 'les_harmonic_mean=.true.,', ''),
                ]
            for label, old, new in guards:
                source = original_source.replace(old, new)
                run_variant(label, args.ranks, " restart_path='split.chk'", expect_success=False)
                assert 'Restart rejected: invalid data or incompatible model/input context' in (work / f'{label}.log').read_text()
                assert not (work / f'{label}.history.csv').exists()
            source = original_source.replace('les_test_ratio=4', 'les_test_ratio=3')
            run_variant('bad_filter', args.ranks, '', expect_success=False)
            assert 'Invalid dynamic LES filter' in (work / 'bad_filter.log').read_text()
            source = original_source
            run_variant('bad_faces', args.ranks, " boundary_faces=6*'outflow', wall_velocity=18*0., restart_path='split.chk'", expect_success=False)
            diagnostic = 'Active Mach number requires inflow face' if args.mach_inlet else 'Restart rejected'
            assert diagnostic in (work / 'bad_faces.log').read_text()
        print(args.les_model, 'PASS: coupled MPI restart, changed ranks, model guards; artifacts:', work)
        return
    schedule = ' injection_events=3\n injection_interval=7.5e-9'
    run_variant('scheduled', args.ranks, schedule)
    schedule_history = rows(work / 'scheduled.history.csv')
    assert [int(row['events']) for row in schedule_history] == [1, 2, 3, 3]
    assert [int(row['parcels']) for row in schedule_history] == [4, 8, 12, 12]
    assert math.isclose(float(schedule_history[-1]['injected_mass']), 3e-10, rel_tol=1e-14)
    assert int(schedule_history[-1]['last_id']) == 12
    assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
               abs(float(row['relative_energy_error'])) < 2e-11 for row in schedule_history)
    run_variant('scheduled_split', args.ranks, schedule + "\n checkpoint_path='scheduled.chk'", stop='7.5e-9')
    split_history = rows(work / 'scheduled_split.history.csv')
    assert int(split_history[-1]['events']) == 2
    run_variant('scheduled_resume', 6 - args.ranks, schedule + "\n restart_path='scheduled.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('scheduled', suffix, keys)
        actual = snapshots('scheduled_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    resumed_history = rows(work / 'scheduled_resume.history.csv')
    assert int(resumed_history[0]['events']) == 2 and int(resumed_history[-1]['events']) == 3
    assert int(resumed_history[-1]['last_id']) == 12 and int(resumed_history[-1]['parcels']) == 12
    profile = (schedule + '\n injection_profile_points=4' +
               '\n injection_profile_times=0,7.5e-9,1.5e-8,2.25e-8' +
               '\n injection_profile_rates=0,0,.02,0')
    run_variant('profiled', args.ranks, profile)
    profile_history = rows(work / 'profiled.history.csv')
    assert [int(row['parcels']) for row in profile_history] == [0, 4, 8, 8]
    assert int(profile_history[-1]['events']) == 3 and int(profile_history[-1]['last_id']) == 8
    assert math.isclose(float(profile_history[-1]['injected_mass']), 1.5e-10, rel_tol=1e-14)
    assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
               abs(float(row['relative_energy_error'])) < 2e-11 for row in profile_history)
    run_variant('profiled_split', args.ranks, profile + "\n checkpoint_path='profiled.chk'", stop='7.5e-9')
    run_variant('profiled_resume', 6 - args.ranks, profile + "\n restart_path='profiled.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('profiled', suffix, keys)
        actual = snapshots('profiled_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    run_variant('profiled_mismatch', args.ranks,
                profile.replace('rates=0,0,.02,0', 'rates=0,0,.03,0') + "\n restart_path='profiled.chk'", expect_success=False)
    assert not (work / 'profiled_mismatch.history.csv').exists()
    ports = """ port_count=2
 port_events=2,2
 port_start=0,7.5e-9
 port_interval=7.5e-9,7.5e-9
 port_number=2,3
 port_mass=1e-10,2e-10
 port_center(:,1)=.002,.002,.002
 port_center(:,2)=.006,.002,.002
 port_axis(:,2)=-1,0,0
 port_fractions(:,1)=1,0
 port_fractions(:,2)=0,1
 port_diameter=2e-5,3e-5
 port_temperature=330,340"""
    run_variant('ports', args.ranks, ports)
    port_history = rows(work / 'ports.history.csv')
    assert [int(row['events']) for row in port_history] == [1, 3, 4, 4]
    assert [int(row['parcels']) for row in port_history] == [2, 7, 10, 10]
    assert math.isclose(float(port_history[-1]['injected_mass']), 6e-10, rel_tol=1e-14)
    assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
               abs(float(row['relative_energy_error'])) < 2e-11 for row in port_history)
    port_particles = snapshots('ports', 'particles', ['id'])
    assert set(port_particles) == {(i,) for i in range(1, 11)}
    for (ident,), row in port_particles.items():
        # At the tie, port 1 reserves IDs 3/4 before port 2 reserves 5/6/7.
        assert float(row['component_mass_2' if ident <= 4 else 'component_mass_1']) == 0
        assert float(row['x']) < .004 if ident <= 4 else float(row['x']) > .004
        assert float(row['u']) > 0 if ident <= 4 else float(row['u']) < 0
    run_variant('ports_split', args.ranks, ports + "\n checkpoint_path='ports.chk'", stop='7.5e-9')
    assert int(rows(work / 'ports_split.history.csv')[-1]['events']) == 3
    run_variant('ports_resume', 6 - args.ranks, ports + "\n restart_path='ports.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('ports', suffix, keys)
        actual = snapshots('ports_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    run_variant('ports_mismatch', args.ranks,
                ports.replace('port_mass=1e-10,2e-10', 'port_mass=1e-10,3e-10') +
                "\n restart_path='ports.chk'", expect_success=False)
    assert not (work / 'ports_mismatch.history.csv').exists()
    run_variant('ports_profile', args.ranks, ports +
                '\n injection_profile_points=4' +
                '\n injection_profile_times=0,7.5e-9,1.5e-8,2.25e-8' +
                '\n injection_profile_rates=0,0,.02,0')
    profiled_ports = rows(work / 'ports_profile.history.csv')
    assert [int(row['parcels']) for row in profiled_ports] == [0, 5, 8, 8]
    assert math.isclose(float(profiled_ports[-1]['injected_mass']), 2.25e-10, rel_tol=1e-14)
    assert int(profiled_ports[-1]['last_id']) == 8
    independent = (ports +
                   '\n injection_profile_points=2' +
                   '\n injection_profile_times=0,1.5e-8' +
                   '\n injection_profile_rates=.01,.01' +
                   '\n port_profile_points=0,3' +
                   '\n port_profile_times(:,2)=7.5e-9,1.5e-8,2.25e-8' +
                   '\n port_profile_rates(:,2)=0,0,.04')
    run_variant('independent', args.ranks, independent)
    independent_history = rows(work / 'independent.history.csv')
    assert [int(row['parcels']) for row in independent_history] == [2, 4, 7, 7]
    assert int(independent_history[-1]['events']) == 4
    assert math.isclose(float(independent_history[-1]['injected_mass']), 3e-10, rel_tol=1e-14)
    assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
               abs(float(row['relative_energy_error'])) < 2e-11 for row in independent_history)
    run_variant('independent_split', args.ranks, independent +
                "\n checkpoint_path='independent.chk'", stop='7.5e-9')
    run_variant('independent_resume', 6 - args.ranks, independent + "\n restart_path='independent.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('independent', suffix, keys)
        actual = snapshots('independent_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    run_variant('independent_mismatch', args.ranks,
                independent.replace('rates(:,2)=0,0,.04', 'rates(:,2)=0,0,.03') +
                "\n restart_path='independent.chk'", expect_success=False)
    assert not (work / 'independent_mismatch.history.csv').exists()
    run_variant('independent_uncovered', args.ranks,
                independent.replace('2.25e-8', '2e-8'), expect_success=False)
    assert not (work / 'independent_uncovered.history.csv').exists()
    # The maximum profile table must survive serialization and checkpoint context.
    maximum = ' port_count=16\n port_events=16*1\n port_profile_points=16*64'
    knots = ','.join(f'{i}e-9' for i in range(64))
    for port in range(1, 17):
        maximum += f'\n port_profile_times(:,{port})={knots}'
    run_variant('maximum_profiles', args.ranks, maximum + "\n checkpoint_path='maximum.chk'", stop='0')
    assert max(map(len, (work / 'maximum.chk').read_text().splitlines())) > 16384
    run_variant('maximum_resume', 6 - args.ranks, maximum + "\n restart_path='maximum.chk'", stop='0')
    maximum_history = rows(work / 'maximum_resume.history.csv')
    assert int(maximum_history[-1]['events']) == 16 and int(maximum_history[-1]['last_id']) == 0
    shaped = (ports + """
 port_events=1,1
 port_start=0,0
 port_nozzle_radius=4e-4,3e-4
 port_nozzle(1)%profile='solid'
 port_nozzle(2)%profile='hollow'
 port_nozzle(2)%swirl=.3
 port_nozzle(2)%spread=.02
 port_size(1)%distribution='lognormal'
 port_size(1)%std=4e-6
 port_size(1)%minimum=1e-5
 port_size(1)%maximum=4e-5
 port_size(2)%shape=2
 port_size(2)%weighting='number'
 port_size(2)%minimum=1e-5
 port_size(2)%maximum=5e-5""")
    run_variant('shaped', args.ranks, shaped)
    run_variant('shaped_initial', args.ranks, shaped + "\n checkpoint_path='shaped.chk'", stop='0')
    shaped_particles = snapshots('shaped_initial', 'particles', ['id'])
    diameters = {1: [], 2: []}
    for (ident,), row in shaped_particles.items():
        port = 1 if ident <= 2 else 2
        rho = 1000 if port == 1 else 800
        d = (6 * float(row['mass']) / (math.pi * rho)) ** (1 / 3)
        assert 1e-5 <= d <= (4e-5 if port == 1 else 5e-5)
        diameters[port].append(d)
        radial = math.hypot(float(row['y']) - .002, float(row['z']) - .002)
        assert 0 < radial <= (4e-4 if port == 1 else 3e-4) * (1 + 1e-12)
        assert math.isclose(float(row['x']), .002 if port == 1 else .006, abs_tol=1e-15)
    assert all(max(ds) > min(ds) * 1.01 for ds in diameters.values())
    run_variant('shaped_resume', 6 - args.ranks, shaped + "\n restart_path='shaped.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('shaped', suffix, keys)
        actual = snapshots('shaped_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    run_variant('shaped_mismatch', args.ranks,
                shaped.replace('swirl=.3', 'swirl=.4') + "\n restart_path='shaped.chk'", expect_success=False)
    assert not (work / 'shaped_mismatch.history.csv').exists()
    run_variant('scalar_shaped', args.ranks, """ nozzle_radius=4e-4
 nozzle%profile='solid'
 size_options%distribution='normal'
 size_options%std=4e-6
 size_options%minimum=1e-5
 size_options%maximum=4e-5""", stop='0')
    scalar_drops = snapshots('scalar_shaped', 'particles', ['id'])
    assert len(scalar_drops) == 4
    assert len({row['mass'] for row in scalar_drops.values()}) == 4
    assert all(math.hypot(float(row['y']) - .002, float(row['z']) - .002) > 0 for row in scalar_drops.values())
    for model in ['etab', 'khrt']:
        breakup = f" breakup_model='{model}'\n surface_tension=.072\n liquid_viscosity=.001\n speed=1000"
        run_variant(model, args.ranks, breakup)
        run_variant(model + '_split', args.ranks, breakup + f"\n checkpoint_path='{model}.chk'", stop='1e-8')
        run_variant(model + '_resume', 6 - args.ranks, breakup + f"\n restart_path='{model}.chk'")
        for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
            reference = snapshots(model, suffix, keys)
            actual = snapshots(model + '_resume', suffix, keys)
            assert actual.keys() == reference.keys()
            for key, expected in reference.items():
                for column, value in expected.items():
                    assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (model, key, column)
        drops = snapshots(model, 'particles', ['id'])
        if model == 'etab':
            assert any(float(row['deformation']) > 0 for row in drops.values())
        else:
            # Initial multiplicity from harmonic mixture density and monodisperse size.
            initial_drop_mass = (1 / (.5 / 1000 + .5 / 800)) * math.pi / 6 * (2e-5) ** 3
            assert all(float(row['multiplicity']) > 1e-10 / (4 * initial_drop_mass) for row in drops.values())
        breakup_history = rows(work / f'{model}.history.csv')
        assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
                   abs(float(row['relative_energy_error'])) < 2e-11 for row in breakup_history)
    run_variant('breakup_mismatch', args.ranks,
                " breakup_model='khrt'\n surface_tension=.073\n liquid_viscosity=.001\n speed=1000\n restart_path='khrt.chk'",
                expect_success=False)
    assert not (work / 'breakup_mismatch.history.csv').exists()
    shedding = (" breakup_model='khrt'\n surface_tension=.072\n liquid_viscosity=.001\n speed=1000" +
                "\n khrt_children=.true.\n khrt_reference_diameter_cubed=1e-25")
    run_variant('shedding', args.ranks, shedding)
    run_variant('shedding_split', args.ranks, shedding + "\n checkpoint_path='shedding.chk'", stop='1e-8')
    run_variant('shedding_resume', 6 - args.ranks, shedding + "\n restart_path='shedding.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('shedding', suffix, keys)
        actual = snapshots('shedding_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    child_drops = snapshots('shedding', 'particles', ['id'])
    assert len(child_drops) > 4 and any(int(row['parent_id']) > 0 for row in child_drops.values())
    assert max(key[0] for key in child_drops) > 4
    child_history = rows(work / 'shedding.history.csv')
    assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
               abs(float(row['relative_energy_error'])) < 2e-11 for row in child_history)
    run_variant('shedding_mismatch', args.ranks,
                shedding.replace('cubed=1e-25', 'cubed=1e-24') + "\n restart_path='shedding.chk'", expect_success=False)
    assert not (work / 'shedding_mismatch.history.csv').exists()
    sgs = ' les_cs=.2\n les_prt=.7\n les_sct=.7\n speed=1000'
    run_variant('sgs_baseline', args.ranks, ' speed=1000')
    run_variant('sgs', args.ranks, sgs)
    run_variant('sgs_split', args.ranks, sgs + "\n checkpoint_path='sgs.chk'", stop='1e-8')
    run_variant('sgs_resume', 6 - args.ranks, sgs + "\n restart_path='sgs.chk'")
    for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
        reference = snapshots('sgs', suffix, keys)
        actual = snapshots('sgs_resume', suffix, keys)
        assert actual.keys() == reference.keys()
        for key, expected in reference.items():
            for column, value in expected.items():
                assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25), (key, column)
    sgs_gas = snapshots('sgs', 'gas', ['i', 'j', 'k'])
    sgs_baseline = snapshots('sgs_baseline', 'gas', ['i', 'j', 'k'])
    assert max(abs(float(sgs_gas[key]['mx']) - float(sgs_baseline[key]['mx']))
               for key in sgs_baseline) > 1e-14
    sgs_history = rows(work / 'sgs.history.csv')
    assert all(abs(float(row['relative_mass_error'])) < 2e-11 and
               abs(float(row['relative_energy_error'])) < 2e-11 for row in sgs_history)
    run_variant('sgs_mismatch', args.ranks,
                sgs.replace('cs=.2', 'cs=.3') + "\n restart_path='sgs.chk'", expect_success=False)
    assert not (work / 'sgs_mismatch.history.csv').exists()
    physical = """
 boundary_faces='outflow','outflow','no_slip_wall','isothermal_wall','periodic','periodic'
 wall_temperature=600,600,600,600,600,600
 wall_velocity(1,3)=2
 center=1e-8,.002,.002
 axis=-1,0,0
 gas_velocity=10,0,0
 les_cs=.2
"""
    run_variant('physical', args.ranks, physical)
    run_variant('physical_split', args.ranks, physical + " checkpoint_path='physical.chk'", stop='1e-8')
    physical_history = rows(work / 'physical.history.csv')
    assert int(physical_history[-1]['parcels']) == 0
    assert abs(float(physical_history[-1]['boundary_energy'])) > 1e-14
    for name, count in [('physical_resume', args.ranks), ('physical_redistributed', 6 - args.ranks)]:
        run_variant(name, count, physical + " restart_path='physical.chk'")
        for suffix, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
            expected = snapshots('physical', suffix, keys)
            actual = snapshots(name, suffix, keys)
            assert expected.keys() == actual.keys()
            for key, record in expected.items():
                for column, value in record.items():
                    if count == args.ranks:
                        assert actual[key][column] == value, (name, key, column)
                    else:
                        assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25)
        resumed_history = rows(work / f'{name}.history.csv')
        for row in resumed_history + physical_history:
            assert abs(float(row['relative_mass_error'])) < 2e-11
            assert abs(float(row['relative_energy_error'])) < 2e-11
        for column in ['boundary_mass', 'boundary_energy']:
            assert math.isclose(float(resumed_history[-1][column]), float(physical_history[-1][column]),
                                rel_tol=2e-11, abs_tol=1e-25)
    run_variant('physical_mismatch', args.ranks,
                physical.replace('600,600,600,600,600,600', '600,600,600,550,600,600') +
                " restart_path='physical.chk'", expect_success=False)
    assert not (work / 'physical_mismatch.history.csv').exists()
    wall_model = """
 boundary_faces='periodic','periodic','periodic','periodic','slip_wall','slip_wall'
 center=.004,.002,1e-10
 axis=0,0,-1
 speed=1
 surface_tension=.072
 liquid_viscosity=.001
 impact%enabled(5)=.true.
 impact%temperature=350,350,350,350,350,350
 impact%evolve=.true.
 impact%liquid_conductivity=.6
 impact%minimum_height=1e-10
"""
    for model_name, model in [('film', wall_model), ('splash', wall_model.replace('speed=1', 'speed=100'))]:
        run_variant(model_name, args.ranks, model)
        run_variant(model_name + '_split', args.ranks, model + f" checkpoint_path='{model_name}.chk'", stop='1e-8')
        reference_history = rows(work / f'{model_name}.history.csv')
        assert int(reference_history[-1]['films']) > 0
        assert float(reference_history[-1]['film_mass']) > 0
        if model_name == 'film':
            assert int(reference_history[-1]['parcels']) == 0
        else:
            assert int(reference_history[-1]['parcels']) > 4
        def film_records(prefix):
            records = [record for path in work.glob(f'{prefix}.rank*.films.csv') for record in rows(path)]
            return sorted(records, key=lambda record: tuple(float(record[key]) for key in
                          ['x', 'y', 'z', 'material_id', 'mass', 'temperature']))
        reference_films = film_records(model_name)
        assert any(float(record['temperature']) > 330 for record in reference_films)
        for record in reference_films:
            assert math.isclose(float(record['mass']), float(record['component_mass_1']) +
                                float(record['component_mass_2']), rel_tol=2e-12)
        for suffix, count in [('resume', args.ranks), ('redistributed', 6 - args.ranks)]:
            name = model_name + '_' + suffix
            run_variant(name, count, model + f" restart_path='{model_name}.chk'")
            for kind, keys in [('gas', ['i', 'j', 'k']), ('particles', ['id'])]:
                expected, actual = snapshots(model_name, kind, keys), snapshots(name, kind, keys)
                assert expected.keys() == actual.keys()
                for key, record in expected.items():
                    for column, value in record.items():
                        if count == args.ranks:
                            assert actual[key][column] == value, (name, key, column)
                        else:
                            assert math.isclose(float(actual[key][column]), float(value), rel_tol=2e-11, abs_tol=1e-25)
            actual_films = film_records(name)
            assert len(actual_films) == len(reference_films)
            for expected, actual in zip(reference_films, actual_films):
                for column, value in expected.items():
                    if count == args.ranks:
                        assert actual[column] == value, (name, column)
                    else:
                        assert math.isclose(float(actual[column]), float(value), rel_tol=2e-11, abs_tol=1e-25)
            for record in rows(work / f'{name}.history.csv') + reference_history:
                assert abs(float(record['relative_mass_error'])) < 2e-11
                assert abs(float(record['relative_energy_error'])) < 2e-11
    run_variant('film_mismatch', args.ranks, wall_model.replace('conductivity=.6', 'conductivity=.7') +
                " restart_path='film.chk'", expect_success=False)
    assert not (work / 'film_mismatch.history.csv').exists()
    print(f'PASS MPI spray application, SGS and restart: {args.ranks} ranks; artifacts: {work}')



if __name__ == '__main__':
    main()
