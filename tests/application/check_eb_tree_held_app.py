#!/usr/bin/env python3
"""Bounded UEINT app lifecycle: regrid, owner migration and changed-rank restart."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial', required=True, type=Path)
    parser.add_argument('--mpi', required=True, type=Path)
    parser.add_argument('--mpiexec', required=True, type=Path)
    parser.add_argument('--work', required=True, type=Path)
    parser.add_argument('--selected', action='store_true')
    parser.add_argument('--separate-momentum', action='store_true')
    parser.add_argument('--subcycling', action='store_true')
    parser.add_argument('--godunov', action='store_true')
    parser.add_argument('--iterations', type=int, default=1)
    parser.add_argument('--chemistry', action='store_true')
    parser.add_argument('--state-redist', action='store_true')
    parser.add_argument('--transport-off', action='store_true')
    parser.add_argument('--embedded-temperature', type=float, default=0)
    parser.add_argument('--embedded-noslip', action='store_true')
    parser.add_argument('--les-model', choices=('smagorinsky', 'wale', 'vreman', 'dynamic'))
    parser.add_argument('--rotation', type=float, default=0)
    parser.add_argument('--smoke-only', action='store_true')
    parser.add_argument('--consistent-velocity', action='store_true')
    parser.add_argument('--paired-pressure', action='store_true')
    parser.add_argument('--connected-directions', action='store_true')
    args = parser.parse_args()
    if args.connected_directions:
        args.separate_momentum = True
    if args.paired_pressure:
        args.consistent_velocity = True
    if args.consistent_velocity:
        args.separate_momentum = True
    momentum_key = 'consistent_velocity' if args.consistent_velocity else 'separate_momentum'
    mismatch_message = 'Restart velocity policy mismatch' if args.consistent_velocity else 'Restart momentum policy mismatch'
    if args.paired_pressure:
        momentum_key = 'paired_pressure'
        mismatch_message = 'Restart pressure policy mismatch'
    if args.connected_directions:
        momentum_key = 'connected_directions'
        mismatch_message = 'Restart connected-direction policy mismatch'
    if args.separate_momentum and (args.godunov or args.state_redist or args.les_model):
        parser.error('Separate momentum requires native held MOL without LES or redistribution')
    if args.godunov and args.rotation:
        parser.error('Rotating-frame integration requires MOL')
    if args.iterations < 1:
        parser.error('iterations must be positive')
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if args.chemistry and not args.selected:
        parser.error('This chemistry diagnostic requires selected pure N2')
    # One root plus a partial child exercises paired transfer. Collapse then
    # recreate it, so a saved file is resumed through an actual topology change.
    template = '''&eb_tree
 dual_energy_enabled=.true.,dual_energy_eta2=2,allow_negative_energy=.false.,
 time_integrator='sdc',sdc_iterations=1,mol_mode=.true.,
 hydro_method='plm_ctu',riemann_solver='pelec',boundary_faces=6*'outflow',
 root_shape=4,4,4,child_count=1,parent=1,ratio=2,
 lower(:,1)=1,1,2,upper(:,1)=2,2,3,
 regular_density=.31,cut_density=.35,initial_pressure=60000,
 geometry_type='plane',plane_axis='z',plane_position=.31,
 time_step=1.e-10,final_time=3.e-10,maximum_steps=20,
 automatic_timestep=.true.,external_forcing=20000,-30000,40000,
 chemistry_enabled=.false.,transport_enabled=.true.,
 regrid_count=2,regrid_times=1.e-10,2.e-10,regrid_child_count=0,1,
 regrid_parent(1,2)=1,regrid_ratio(1,2)=2,
 regrid_lower(:,1,2)=1,1,2,regrid_upper(:,1,2)=2,2,3,
 output_file='result.csv',
'''
    template = template.replace('sdc_iterations=1', f'sdc_iterations={args.iterations}')
    if args.godunov:
        template = template.replace('mol_mode=.true.', 'mol_mode=.false.')
    if args.separate_momentum:
        template += ' separate_momentum=.true.,\n'
    if args.consistent_velocity:
        template += ' consistent_velocity=.true.,\n'
    if args.paired_pressure:
        template += '\n paired_pressure=.true.,\n'
    if args.connected_directions:
        template += '\n connected_directions=.true.,\n'
    if args.selected:
        template += " composition_count=1,composition_species='N2',composition_mole_fractions=1,\n"
    if args.chemistry:
        template = template.replace('chemistry_enabled=.false.', 'chemistry_enabled=.true.')
    if args.transport_off:
        template = template.replace('transport_enabled=.true.', 'transport_enabled=.false.')
    if args.subcycling:
        template += ' subcycling=.true.,\n'
    if args.state_redist:
        template += " redistribution='state_redist',state_redist_target_volume_fraction=.9,\n"
    if args.les_model:
        # Fully regular cells and reflecting physical walls produce a nonzero
        # LES strain after the forced predictor. Cut-cell LES remains rejected.
        template = template.replace("plane_position=.31", "plane_position=-100")
        template = template.replace("boundary_faces=6*'outflow'", "boundary_faces=6*'no_slip_wall'")
        template = template.replace('1.e-10', '1.e-7').replace('2.e-10', '2.e-7').replace('3.e-10', '3.e-7')
    if args.embedded_temperature:
        template += f' embedded_wall_temperature={args.embedded_temperature:.17e},\n'
    if args.embedded_noslip:
        template += ' embedded_noslip=.true.,\n'
    if args.rotation:
        template += f' rotating_frame_enabled=.true.,rotation_angular_speed={args.rotation:.17e},\n'
        if abs(args.rotation) >= 1.e8:
            template = template.replace('1.e-10', '1.e-25').replace('2.e-10', '2.e-25').replace('3.e-10', '3.e-25')
    env = os.environ.copy()
    env['OMPI_MCA_rmaps_base_oversubscribe'] = '1'
    results = []

    def run(name, ranks=0, extra='', success=True, source=None, les_settings=None):
        directory = work / name
        directory.mkdir(exist_ok=True)
        text = template if source is None else source
        les_input = '' if not args.les_model else f'&eb_tree_les les_enabled=.true.,les_model=\'{args.les_model}\',/\n'
        if les_settings is not None:
            les_input = les_settings
        (directory / 'input.nml').write_text(text + extra + '\n/\n' + les_input)
        command = [str(args.serial.resolve()), 'input.nml']
        if ranks:
            command = [str(args.mpiexec.resolve()), '-n', str(ranks), str(args.mpi.resolve()), 'input.nml']
        result = subprocess.run(command, cwd=directory, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=240)
        (directory / 'run.log').write_text(result.stdout)
        binary = args.mpi if ranks else args.serial
        results.append({'name': name, 'command': command, 'returncode': result.returncode,
                        'executable_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
                        'input_sha256': hashlib.sha256((directory / 'input.nml').read_bytes()).hexdigest()})
        (work / 'commands.json').write_text(json.dumps(results, indent=2) + '\n')
        if (result.returncode == 0) != success:
            raise AssertionError(f'{directory}: exit {result.returncode}\n{result.stdout}')
        return directory

    serial = run('serial')
    if args.smoke_only:
        print(f'PASS: held serial app smoke {serial}')
        return
    if args.les_model:
        disabled = run('les_disabled', les_settings='')
        if args.les_model == 'smagorinsky':
            assert (disabled / 'result.csv').read_bytes() != (serial / 'result.csv').read_bytes(), \
                'LES must change the forced reflecting-wall state'
    full = run('mpi2', 2, 'rebalance_interval_steps=1,rebalance_minimum_gain=0.,')
    third = run('mpi3', 3)
    stop = "checkpoint_file='saved.chk',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,"
    stopped = run('stopped2', 2, stop)
    assert 'Regrid stage: 1 nodes: 1' in (stopped / 'run.log').read_text()
    if args.separate_momentum:
        checkpoint = stopped / 'saved.chk'
        token = ('PELEF_EB_TREE_CONNECTED_DIRECTIONS_V1' if args.connected_directions else
                         'PELEF_EB_TREE_PAIRED_PRESSURE_V1' if args.paired_pressure else
                         'PELEF_EB_TREE_CONSISTENT_VELOCITY_V1' if args.consistent_velocity else 'PELEF_EB_TREE_SEPARATE_MOMENTUM_V1')
        assert checkpoint.read_text().startswith(token + '\n')
        old_template = template.replace(momentum_key + '=.true.', momentum_key + '=.false.')
        old = run('legacy_stopped2', 2, stop, source=old_template)
        assert token not in (old / 'saved.chk').read_text()
        for name, source, saved in [('candidate_to_legacy', old_template, checkpoint),
                                    ('legacy_to_candidate', template, old / 'saved.chk')]:
            before = saved.read_bytes()
            denied = run(name, 3, f"restart_file='{saved}',", success=False, source=source)
            assert mismatch_message in (denied / 'run.log').read_text()
            assert not list(denied.glob('result.csv*'))
            assert saved.read_bytes() == before
        for name, source in [('missing_held', template.replace('dual_energy_enabled=.true.',
                                                               'dual_energy_enabled=.false.')),
                             ('not_mol', template.replace('mol_mode=.true.', 'mol_mode=.false.')),
                             ('prepared_redist', template + " redistribution='state_redist',")]:
            for ranks in (0, 2):
                denied = run(f'{name}_{ranks}', ranks, success=False, source=source)
                assert 'Separate momentum requires held MOL without redistribution' in (denied / 'run.log').read_text()
                assert not list(denied.glob('result.csv*'))

    if args.connected_directions:
        for ranks in (0, 2):
            denied = run(f'connected_requires_momentum_{ranks}', ranks, success=False,
                         source=template.replace('separate_momentum=.true.', 'separate_momentum=.false.'))
            assert 'Connected directions requires separate momentum' in (denied / 'run.log').read_text()
            assert not list(denied.glob('result.csv*'))

    if args.paired_pressure:
        for name, source, message in [
            ('velocity_required', template.replace('consistent_velocity=.true.', 'consistent_velocity=.false.'),
             'requires consistent velocity'),
            ('rotating_pressure', template + ' rotating_frame_enabled=.true.,', 'requires a nonrotating frame')]:
            for ranks in (0, 2):
                denied = run(f'{name}_{ranks}', ranks, success=False, source=source)
                assert message in (denied / 'run.log').read_text()
                assert not list(denied.glob('result.csv*'))

    restart = f"restart_file='{stopped / 'saved.chk'}',"
    resumed = run('resumed3', 3, restart)
    to_serial = run('mpi_to_serial', 0, restart)
    serial_stop = run('serial_stop', 0, stop)
    from_serial = run('serial_to_mpi', 3, f"restart_file='{serial_stop / 'saved.chk'}',")
    assert 'PELEF_EB_TREE_CHECKPOINT_V10' in (stopped / 'saved.chk').read_text()
    assert 'PELEF_EB_TREE_CHECKPOINT_V9' in (serial_stop / 'saved.chk').read_text()
    with (stopped / 'result.csv.budget.csv').open() as stream:
        stopped_row = list(csv.DictReader(stream))[4]
        stopped_energy = stopped_row['numerical_energy']
        stopped_pressure_work = stopped_row['internal_energy_pressure_work']
    assert math.isfinite(float(stopped_pressure_work))
    if args.godunov:
        assert float(stopped_pressure_work) != 0, 'fixture must exercise actual Godunov UEINT pressure work'
    else:
        assert float(stopped_pressure_work) == 0, 'MOL has no Godunov pdivu source'
    assert float(stopped_energy) != 0, 'fixture must create a real numerical correction before restart'
    stop_time = 1.e-25 if abs(args.rotation) >= 1.e8 else (1.e-7 if args.les_model else 1.e-10)
    read_only = run('restore_only3', 3, restart + f'final_time={stop_time:.17e},')
    with (read_only / 'result.csv.budget.csv').open() as stream:
        restored = list(csv.DictReader(stream))[4]
        assert restored['numerical_energy'] == stopped_energy, 'exact global prefix restore'
        assert restored['internal_energy_pressure_work'] == stopped_pressure_work, 'exact pressure-work prefix restore'
    thermo = json.loads((Path(__file__).resolve().parents[2] / 'mechanisms/h2o2_full.json').read_text())['thermo']
    n2 = next(item for item in thermo if item['name'] == 'N2')
    expected = (serial / 'result.csv').read_bytes()
    worst = 0
    for directory in (serial, full, third, resumed, to_serial, from_serial):
        assert (directory / 'result.csv').read_bytes() == expected, directory
        with (directory / 'result.csv').open() as stream:
            cells = list(csv.DictReader(stream))
        with (directory / 'result.csv.budget.csv').open() as stream:
            budget = list(csv.DictReader(stream))
        assert {int(cell['level']) for cell in cells} == {0, 1}
        assert all(math.isfinite(float(cell['rho_internal_energy_J_m3'])) for cell in cells)
        assert all(math.isfinite(float(cell['temperature'])) and float(cell['temperature']) > 0 for cell in cells)
        if args.selected:
            # Independent NASA polynomial at the reported temperature, using
            # the byte-pinned bundle rather than E minus kinetic/potential work.
            for cell in cells:
                t = float(cell['temperature'])
                a = n2['low_coefficients'] if t < n2['temperature_mid'] else n2['high_coefficients']
                e = 8314.46261815324/n2['molecular_weight'] * (
                    (a[0]-1)*t + a[1]*t*t/2 + a[2]*t**3/3 + a[3]*t**4/4 + a[4]*t**5/5 + a[5])
                actual, reference = float(cell['rho_internal_energy_J_m3']), float(cell['rho'])*e
                assert abs(actual-reference) <= 3.e-10*max(1, abs(actual), abs(reference)), directory
        fields = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE'] + [f for f in cells[0] if f.startswith('rhoY_')]
        assert len(fields) == len(budget)
        for field, row in zip(fields, budget):
            initial, final = float(row['initial']), float(row['final'])
            inventory = math.fsum(float(cell[field]) * float(cell['fluid_volume']) for cell in cells)
            assert abs(inventory-final) <= 3.e-14 * max(1, abs(final)), (directory, field)
            amounts = [float(row[key]) for key in ('boundary', 'wall', 'chemistry', 'forcing', 'numerical_energy')]
            residual = final - initial - math.fsum(amounts)
            scale = max(1, abs(initial), abs(final), math.fsum(abs(x) for x in amounts))
            assert abs(residual)/scale <= 5.e-10, (directory, field, residual/scale)
            assert abs(residual-float(row['residual'])) <= 3.e-14*scale
            if field != 'rhoE':
                assert float(row['numerical_energy']) == 0, 'correction belongs only to total energy'
            worst = max(worst, abs(residual)/scale)
    # A continuation must retain the global numerical prefix, rather than
    # reconstructing it from the final gas/scalar or adding per-patch copies.
    # Compare each continuation with its native destination execution. Serial
    # and MPI may group the composite scalar sum differently by one final bit.
    for name, destination in ((resumed, third), (to_serial, serial), (from_serial, third)):
        with (name / 'result.csv.budget.csv').open() as stream:
            observed = list(csv.DictReader(stream))[4]
        with (destination / 'result.csv.budget.csv').open() as stream:
            reference = list(csv.DictReader(stream))[4]
        assert observed['numerical_energy'] == reference['numerical_energy'], name
    for name, extra in [('wrong_eta2', 'dual_energy_eta2=3,'),
                        ('wrong_small', 'allow_small_energy=.false.,'),
                        ('wrong_negative', 'allow_negative_energy=.true.,'),
                        ('wrong_update', 'update_total_from_internal=.false.,')]:
        directory = work / name
        directory.mkdir(exist_ok=True)
        paths = [directory / f for f in ('result.csv', 'result.csv.budget.csv', 'protected.chk')]
        for path in paths:
            path.write_bytes(b'preserve\n')
        run(name, 3, restart + extra + "checkpoint_file='protected.chk',", False)
        assert all(path.read_bytes() == b'preserve\n' for path in paths)
    legacy = template.replace('allow_negative_energy=.false.,', '').replace(
        'dual_energy_enabled=.true.', 'dual_energy_enabled=.false.').replace('dual_energy_eta2=2',
                                                                                             'dual_energy_eta2=1.e-4')
    legacy = legacy.replace('separate_momentum=.true.', 'separate_momentum=.false.')
    legacy = legacy.replace('consistent_velocity=.true.', 'consistent_velocity=.false.')
    legacy = legacy.replace('paired_pressure=.true.', 'paired_pressure=.false.')
    legacy = legacy.replace('connected_directions=.true.', 'connected_directions=.false.')
    run('held_to_legacy', 3, restart, False, legacy)
    if abs(args.rotation) < 1.e8:
        old_serial = run('legacy_serial', source=legacy)
        old_mpi = run('legacy_mpi3', 3, source=legacy)
        assert (old_serial / 'result.csv').read_bytes() == (old_mpi / 'result.csv').read_bytes()
        assert 'rho_internal_energy_J_m3' not in (old_serial / 'result.csv').read_text().splitlines()[0]
        assert 'numerical_energy' not in (old_serial / 'result.csv.budget.csv').read_text().splitlines()[0]
    # A nonfinite new scalar history refuses before touching protected outputs.
    damaged = (stopped / 'saved.chk').read_text().splitlines(keepends=True)
    header = next(i for i, line in enumerate(damaged) if line.strip() == 'PELEF_EB_TREE_CHECKPOINT_V10')
    nodes = int(damaged[header + 2])
    prefix_line = header + 5 + nodes
    raw = damaged[prefix_line].removesuffix('\n')
    assert len(raw.split()) == 6
    damaged[prefix_line] = raw[:-28] + ' ' * 24 + 'NaN \n'
    assert len(damaged[prefix_line]) == len(raw) + 1
    bad_prefix = work / 'nonfinite_pressure_work.chk'
    bad_prefix.write_text(''.join(damaged))
    directory = work / 'nonfinite_pressure_restart'
    directory.mkdir(exist_ok=True)
    paths = [directory / f for f in ('result.csv', 'result.csv.budget.csv', 'protected.chk')]
    for path in paths:
        path.write_bytes(b'preserve\n')
    run('nonfinite_pressure_restart', 3,
        f"restart_file='{bad_prefix}',checkpoint_file='protected.chk',", False)
    assert all(path.read_bytes() == b'preserve\n' for path in paths)
    print(f'PASS: held app lifecycle, exact serial/MPI/restart gas/T/UEINT and scalar prefix restore; '
          f'inventory residual {worst:.3e}')


if __name__ == '__main__':
    main()
