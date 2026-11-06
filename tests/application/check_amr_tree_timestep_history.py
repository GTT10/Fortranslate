"""Check common/recursive tree timestep policies, regrid hold and rank-changing restart."""
import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile


def header(path):
    records = path.read_text().splitlines()
    if records[0] in ('PELEF_AMR_TREE_INTERNAL_ENERGY_V1', 'PELEF_AMR_TREE_INTERNAL_ENERGY_V2'):
        values = [float(value) for value in records[1].split()]
        assert len(values) == (2 if records[0].endswith('V2') else 1)
        assert all(math.isfinite(value) for value in values)
        records = records[2:]
    assert records[0] == 'PELEF_AMR_TREE_TIMESTEP_HISTORY_V1'
    length = int(records[2]); cursor = 3; seen = 0
    while seen < length:
        seen += len(records[cursor]) + 1; cursor += 1
    assert seen == length
    progress = records[cursor + 1].split()
    history = records[cursor + 2].split()
    return float(progress[0]), int(progress[4]), float(history[0]), history[1] == 'T'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('serial', 'mpi', 'mpiexec', 'work'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--adaptive-only', action='store_true')
    parser.add_argument('--mol', action='store_true')
    parser.add_argument('--wall-transport', action='store_true')
    parser.add_argument('--mol-default', action='store_true')
    parser.add_argument('--dual-energy', action='store_true')
    args = parser.parse_args()
    args.mol = args.mol or args.mol_default
    if args.wall_transport:
        assert args.mol and not args.mol_default
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='amr_timestep_', dir=args.work.resolve()))
    tree = work / 'tree.dat'; tree.write_text('1\n1 1 2 1 4 1 4 2\n')
    coarse = work / 'coarse.dat'; coarse.write_text('0\n')
    schedule = work / 'schedule.dat'; schedule.write_text(f"1\n5e-9 '{coarse}'\n")
    variants = [('auto', False, False), ('auto', True, False), ('auto', False, True),
                ('auto', True, True), ('fixed', False, False), ('combined', False, False), ('adaptive', False, False)]
    if args.mol:
        variants = variants[:4]
    if args.adaptive_only:
        variants = [('adaptive', False, False)]
    for parallel in (False, True):
        for policy, recursive, regrid in variants:
            name = f"{'mpi' if parallel else 'serial'}_{policy}_{int(recursive)}_{int(regrid)}"
            first, last, steps = 5e-9, 6.05e-9, 3
            extra = ''
            if regrid:
                last, steps = 1.05e-9, 4
            if policy == 'adaptive':
                last, steps = 1.55e-9, 4
            if policy == 'fixed':
                first, last, steps = 3e-9, 1.55e-9, 4
                extra = 'fixed_dt=6e-9,maximum_dt=6e-9,'
            elif policy == 'combined':
                first, last, steps = 8e-9, 2.55e-9, 3
                extra = 'initial_dt=8e-9,fixed_dt=6e-9,maximum_dt=6e-9,'
            else:
                extra = 'maximum_dt=1e-8,'
            config = f"""&reactive_3d
nx=4,ny=4,nz=4,x_upper=.004,y_upper=.004,z_upper=.004
final_time=1.655e-8,maximum_steps=100,cfl=.2
problem='uniform_reactor',thermo_model='elementary',x_h2=0,x_o2=0,x_h=0,x_o=0,x_oh=0,x_n2=1
initial_temperature=800,initial_pressure=101325
reconstruction='characteristic_ppm',riemann_solver='pelec',boundary_condition='periodic'
time_integrator='sdc',sdc_iterations=2,chemistry_enabled=.false.,transport_enabled=.false.
output_file='gas.csv' /
&amr_tree_timestepping subcycling={'.true.' if recursive else '.false.'},
{extra}initial_timestep_factor=.5,maximum_timestep_growth=1.1 /
"""

            if args.mol:
                config = config.replace("reconstruction='characteristic_ppm'", "reconstruction='characteristic_plm'")
                config = config.replace("time_integrator='sdc'", "time_integrator='mol',mol_lax_friedrichs=.false.")

            if args.wall_transport:
                config = config.replace("boundary_condition='periodic'",
                                        "boundary_condition='no_slip_wall',initial_velocity_x=1,"
                                        "initial_velocity_y=2,initial_velocity_z=3")
                config = config.replace('transport_enabled=.false.', 'transport_enabled=.true.')

            if args.dual_energy:
                assert args.mol
                config = config.replace("time_integrator='mol'",
                                        "time_integrator='mol',dual_energy_eta2=1,allow_negative_energy=.false.")

            if args.mol_default:
                config = config.replace('sdc_iterations=2,', '')
                config = config.replace("problem='uniform_reactor'", "problem='entropy_wave',density_wave_amplitude=.01,"
                                        "initial_velocity_x=10,initial_velocity_y=20,initial_velocity_z=30")

            if policy == 'adaptive':
                config += '&amr_tree_adaptation regrid_interval_steps=1,refinement_levels=0 /\n'

            def run(label, control='', ranks=1, reject=False, changed=False):
                folder = work / f'{name}_{label}'; folder.mkdir()
                text = config.replace('maximum_timestep_growth=1.1', 'maximum_timestep_growth=1.2') if changed else config
                (folder / 'input.nml').write_text(text + f'&amr_tree_checkpoint {control} /\n')
                command = ([str(args.mpiexec.resolve()), '-n', str(ranks), str(args.mpi.resolve())]
                           if parallel else [str(args.serial.resolve())]) + ['input.nml', str(coarse if policy == 'adaptive' else tree)]
                if regrid:
                    command.append(str(schedule))
                result = subprocess.run(command, cwd=folder, text=True, capture_output=True, timeout=120)
                (folder / 'run.log').write_text(result.stdout + result.stderr)
                if reject:
                    assert result.returncode != 0, result.stdout + result.stderr
                else:
                    assert result.returncode == 0, result.stdout + result.stderr
                return folder

            full = run('full', "checkpoint_file='full.chk',checkpoint_interval_steps=1,checkpoint_indexed=.true.",
                       ranks=2 if parallel else 1)
            stopped = run('stop', "checkpoint_file='stop.chk',checkpoint_interval_steps=1,"
                          "checkpoint_indexed=.true.,stop_after_checkpoint=.true.", ranks=2 if parallel else 1)
            restart = f"restart_file='{stopped / 'stop.chk'}'"
            resumed = run('resume', restart + ",checkpoint_file='resumed.chk',checkpoint_interval_steps=1,"
                          "checkpoint_indexed=.true.", ranks=1)
            time, count, dt, held = header(full / 'full.chk')
            assert count == steps and abs(time-1.655e-8) < 1e-23 and abs(dt-last) < 1e-23, (name, time, count, dt)
            time, count, dt, held = header(stopped / 'stop.chk')
            assert count == 1 and time == first and dt == first and held == (regrid or policy == 'adaptive'), (name, time, count, dt, held)
            assert header(full / 'full.chk') == header(resumed / 'resumed.chk'), name + ' history'
            assert (full / 'gas.csv').read_bytes() == (resumed / 'gas.csv').read_bytes(), name + ' gas output'
            assert (full / 'full.chk').read_bytes() == (resumed / 'resumed.chk').read_bytes(), name + ' snapshot'
            if args.mol_default:
                with (full / 'gas.csv.budget.csv').open() as stream:
                    for row in csv.DictReader(stream):
                        initial, final, boundary, source = (float(row[key]) for key in
                                                           ('initial', 'final', 'boundary_change', 'source_change'))
                        repair = float(row.get('energy_correction', 0))
                        assert abs(final-initial-boundary-source-repair) <= 2e-10*max(1e-12, abs(initial), abs(final)), name
                # The omitted value must run the same native one-pass algorithm,
                # including serialization, as an explicit one-iteration request.
                config = config.replace("time_integrator='mol'", "time_integrator='mol',sdc_iterations=1")
                explicit = run('explicit_one', "checkpoint_file='explicit.chk',checkpoint_interval_steps=1,"
                               "checkpoint_indexed=.true.", ranks=2 if parallel else 1)
                assert (full / 'gas.csv').read_bytes() == (explicit / 'gas.csv').read_bytes(), name + ' default gas'
                assert (full / 'full.chk').read_bytes() == (explicit / 'explicit.chk').read_bytes(), name + ' default snapshot'

            run('changed_growth', restart, ranks=2 if parallel else 1, reject=True, changed=True)
            if args.dual_energy:
                config = config.replace('dual_energy_eta2=1,', 'dual_energy_eta2=.5,')
                run('changed_energy', restart, ranks=2 if parallel else 1, reject=True)
                config = config.replace('dual_energy_eta2=.5,', 'dual_energy_eta2=1,')
            if args.mol:
                config = config.replace('mol_lax_friedrichs=.false.', 'mol_lax_friedrichs=.true.')
                run('changed_lf', restart, ranks=2 if parallel else 1, reject=True)
            if policy == 'combined':
                records = (stopped / 'stop.chk').read_text().splitlines()
                offset = 2 if records[0].startswith('PELEF_AMR_TREE_INTERNAL_ENERGY_V') else 0
                assert records[offset] == 'PELEF_AMR_TREE_TIMESTEP_HISTORY_V1'
                cursor, seen = offset + 3, 0
                while seen < int(records[offset + 2]):
                    seen += len(records[cursor]) + 1; cursor += 1
                legacy = work / f'{name}_missing_history.chk'
                legacy.write_text('\n'.join(line for i, line in enumerate(records)
                                           if i not in (offset, cursor + 2)) + '\n')
                run('missing_history', f"restart_file='{legacy}'", ranks=2 if parallel else 1, reject=True)

    print('PASS', 'MOL default' if args.mol_default else ('MOL' if args.mol else 'SDC'),
          'tree timestep history, common/recursive, post-regrid hold and MPI2-to1 restart:', work)


if __name__ == '__main__':
    main()
