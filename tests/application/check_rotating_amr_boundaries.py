"""Rotating MOL tree physical boundaries, budgets and changed-rank restart.

This is an implementation check on a small stationary regular mesh. It does
not establish fixed-PeleC trajectory parity or physical validation.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('serial', 'mpi', 'mpiexec', 'work'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='rotating_boundaries_', dir=args.work.resolve()))
    print(work, flush=True)
    tree = work / 'tree.dat'
    tree.write_text('2\n1 1 1 1 1 1 1 2\n2 1 1 1 1 1 1 2\n')
    report = {'scope': __doc__, 'executables': {
        key: {'path': str(getattr(args, key).resolve()),
              'sha256': hashlib.sha256(getattr(args, key).read_bytes()).hexdigest()}
        for key in ('serial', 'mpi')}, 'cases': []}
    for boundary in ('slip_wall', 'no_slip_wall', 'isothermal_wall', 'pressure'):
        for recursive in (False, True):
            name = f'{boundary}_sub{int(recursive)}'
            faces = [boundary, boundary, 'outflow', 'outflow', 'periodic', 'periodic']
            controls = ''
            if boundary == 'pressure':
                faces[:2] = ['outflow', 'outflow']
                controls = 'outlet_pressure=0,101000,0,0,0,0,'
            if boundary == 'no_slip_wall':
                controls = 'wall_velocity=0,3,0,0,3,0,12*0,'
            config = f"""&reactive_3d
time_integrator='mol',sdc_iterations=2,dual_energy_eta2=1,allow_negative_energy=.false.,
nx=4,ny=4,nz=4,x_lower=-.002,x_upper=.002,y_lower=.003,y_upper=.009,z_upper=.008,
final_time=3e-10,cfl=.1,maximum_steps=100,problem='uniform_reactor',thermo_model='full_h2o2',
initial_temperature=900,initial_pressure=100000,
initial_velocity_x=10,initial_velocity_y=20,initial_velocity_z=-5,
reconstruction='characteristic_plm',riemann_solver='pelec',boundary_condition='outflow',
boundary_faces={','.join(repr(face) for face in faces)},wall_temperature=6*850,{controls}
chemistry_enabled=.false.,transport_enabled=.true.,output_file='gas.csv' /
&reactive_rotation enabled=.true.,axis=2,angular_speed=400,axis_point=.0004,-.0007,0 /
&amr_tree_timestepping subcycling={'.true.' if recursive else '.false.'},maximum_dt=1e-10 /
"""

            def run(label, ranks, restart='', stop=False, changed=False):
                directory = work / f'{name}_{label}'
                directory.mkdir()
                checkpoint = directory / 'state.chk'
                text = config.replace('angular_speed=400', 'angular_speed=401') if changed else config
                text += f"&amr_tree_checkpoint checkpoint_file='{checkpoint}',checkpoint_interval_steps=1," \
                        f"checkpoint_indexed=.true.,restart_file='{restart}'," \
                        f"stop_after_checkpoint={'.true.' if stop else '.false.'} /\n"
                (directory / 'input.nml').write_text(text)
                command = ([str(args.mpiexec.resolve()), '-n', str(ranks), str(args.mpi.resolve())]
                           if ranks else [str(args.serial.resolve())]) + ['input.nml', str(tree)]
                result = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=120)
                log = result.stdout + result.stderr
                (directory / 'run.log').write_text(log)
                if changed:
                    assert result.returncode != 0 and 'restart' in log.lower(), log
                    assert not (directory / 'gas.csv').exists()
                    return directory
                assert result.returncode == 0, log
                if stop:
                    assert 'Stopped after checkpoint at step 1' in log, log
                else:
                    assert 'Steps: 3' in log and 'Final time:' in log, log
                return directory

            serial = run('serial', 0)
            parallel = run('mpi', 2)
            stopped = run('stop', 2, stop=True)
            resumed = run('resume', 3, restart=stopped / 'state.chk')
            assert (serial / 'gas.csv').read_bytes() == (parallel / 'gas.csv').read_bytes(), name + ' serial/MPI state'
            assert (parallel / 'gas.csv').read_bytes() == (resumed / 'gas.csv').read_bytes(), name + ' restart state'
            assert (parallel / 'state.chk').read_bytes() == (resumed / 'state.chk').read_bytes(), name + ' restart snapshot'
            maximum = 0.
            for directory in (serial, parallel, resumed):
                with (directory / 'gas.csv.budget.csv').open() as stream:
                    for row in csv.DictReader(stream):
                        initial, final, b, source = (float(row[key]) for key in
                                                   ('initial', 'final', 'boundary_change', 'source_change'))
                        correction = float(row.get('energy_correction', 0))
                        residual = abs(final-initial-b-source-correction) / max(1e-12, abs(initial), abs(final))
                        assert residual <= 3e-12, (name, row, residual)
                        maximum = max(maximum, residual)
            run('changed_frame', 3, restart=stopped / 'state.chk', changed=True)
            report['cases'].append({'name': name, 'maximum_budget_residual': maximum})
            (work / 'checks.json').write_text(json.dumps(report, indent=2) + '\n')
            print(name, maximum, flush=True)


if __name__ == '__main__':
    main()
