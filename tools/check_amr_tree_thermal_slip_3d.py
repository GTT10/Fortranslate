#!/usr/bin/env python3
"""Check explicit tree thermal-slip controls, MPI identity and restart binding."""
import argparse
import csv
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('serial', 'mpi', 'work'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--mpiexec', default='/usr/bin/mpiexec')
    parser.add_argument('--selected', action='store_true')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='thermal_slip_', dir=args.work.resolve()))
    tree = work / 'tree.dat'
    tree.write_text('1\n1 1 2 1 4 1 4 2\n')
    env = dict(os.environ, OMP_NUM_THREADS='1', OMPI_MCA_rmaps_base_oversubscribe='1')
    runs = 0
    for integrator in ('sdc', 'mol'):
        for recursive in (False, True):
            case = work / f'{integrator}_{int(recursive)}'
            case.mkdir()
            config = f"""&reactive_3d
nx=4,ny=4,nz=4,x_upper=.004,y_upper=.004,z_upper=.004
final_time=1.655e-8,maximum_steps=100,cfl=.2
problem='uniform_reactor',thermo_model='elementary',x_h2=0,x_o2=0,x_h=0,x_o=0,x_oh=0,x_n2=1
initial_temperature=800,initial_pressure=101325
reconstruction='characteristic_plm',riemann_solver='pelec',boundary_condition='slip_wall'
initial_velocity_x=0,initial_velocity_y=2,initial_velocity_z=3
wall_temperature=500
time_integrator='{integrator}',chemistry_enabled=.false.,transport_enabled=.true.
output_file='gas.csv' /
&amr_tree_timestepping subcycling={'.true.' if recursive else '.false.'},
maximum_dt=1e-8,initial_timestep_factor=.5,maximum_timestep_growth=1.1,
isothermal_slip_faces=.true.,.true.,.false.,.false.,.false.,.false. /
"""
            if integrator == 'sdc':
                config = config.replace("reconstruction='characteristic_plm'", "reconstruction='characteristic_ppm'")
                config = config.replace("time_integrator='sdc'", "time_integrator='sdc',sdc_iterations=2")
            if args.selected:
                config = config.replace("thermo_model='elementary'", "thermo_model='selected',"
                                        "composition_count=1,composition_species='N2',composition_mole_fractions=1")

            def run(label, content, ranks=0, control='', error=None):
                nonlocal runs
                folder = case / label
                folder.mkdir()
                (folder / 'input.nml').write_text(content + f'&amr_tree_checkpoint {control} /\n')
                command = ([args.mpiexec, '-n', str(ranks), str(args.mpi.resolve())] if ranks
                           else [str(args.serial.resolve())]) + ['input.nml', str(tree)]
                result = subprocess.run(command, cwd=folder, env=env, capture_output=True, text=True, timeout=180)
                log = result.stdout + result.stderr
                (folder / 'run.log').write_text(log)
                runs += 1
                if error:
                    assert result.returncode != 0 and error in log, (label, log)
                    assert not (folder / 'gas.csv').exists(), label
                else:
                    assert result.returncode == 0, (label, log)
                return folder

            full = run('serial_full', config)
            parallel = run('mpi_full', config, 2)
            outputs = ('gas.csv', 'gas.csv.budget.csv')

            def same(folder):
                for name in outputs:
                    assert (folder / name).read_bytes() == (full / name).read_bytes(), (folder, name)

            same(parallel)
            with (full / 'gas.csv.budget.csv').open() as stream:
                for row in csv.DictReader(stream):
                    initial, final, boundary, source = (float(row[key]) for key in
                                                       ('initial', 'final', 'boundary_change', 'source_change'))
                    repair = float(row.get('energy_correction', 0))
                    assert abs(final-initial-boundary-source-repair) <= 2e-10*max(1e-12, abs(initial), abs(final)), row
            save = "checkpoint_file='stop.chk',checkpoint_interval_steps=1,stop_after_checkpoint=.true."
            stopped = run('mpi_stop', config, 2, save)
            serial_stop = run('serial_stop', config, 0, save)
            assert (stopped / 'stop.chk').read_bytes() == (serial_stop / 'stop.chk').read_bytes()
            restart = f"restart_file='{stopped / 'stop.chk'}'"
            for ranks in (0, 1, 3):
                same(run(f'resume_{ranks}', config, ranks, restart))
            for ranks in (0, 3):
                changed = config.replace('isothermal_slip_faces=.true.,.true.', 'isothermal_slip_faces=.false.,.true.')
                run(f'changed_mask_{ranks}', changed, ranks, restart, 'restart checkpoint/context')
                changed = config.replace('wall_temperature=500', 'wall_temperature=501')
                run(f'changed_temperature_{ranks}', changed, ranks, restart, 'restart checkpoint/context')
            adiabatic = config.replace('isothermal_slip_faces=.true.,.true.', 'isothermal_slip_faces=.false.,.false.')
            cold = run('adiabatic_explicit', adiabatic)
            omitted = '\n'.join(line for line in adiabatic.splitlines() if not line.startswith('isothermal_slip_faces=')) + '\n/\n'
            default = run('adiabatic_default', omitted)
            assert (cold / 'gas.csv').read_bytes() == (default / 'gas.csv').read_bytes()
            assert (cold / 'gas.csv').read_bytes() != (full / 'gas.csv').read_bytes(), 'Thermal mask has no effect'
            for ranks in (0, 3):
                run(f'wrong_face_{ranks}', config.replace("boundary_condition='slip_wall'", "boundary_condition='outflow'"),
                    ranks, error='Thermal slip requires slip_wall faces')
                run(f'no_conduction_{ranks}', config.replace('transport_enabled=.true.', 'transport_enabled=.false.'),
                    ranks, error='Thermal slip requires SDC/MOL conduction')
    print(f'PASS: {runs} thermal-slip application runs, common/recursive SDC/MOL, MPI identity, restart and refusal: {work}')


if __name__ == '__main__':
    main()
