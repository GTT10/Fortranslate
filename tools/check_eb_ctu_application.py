#!/usr/bin/env python3
"""Check single-level EB CTU composition and method-bound restart."""
import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--selected', action='store_true', help='Use the selected H2O2 fixture')
    parser.add_argument('--inflow', action='store_true', help='Exercise prescribed physical inflow and restart binding')
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='run-', dir=args.work_dir.resolve()))
    base = """&reactive_eb_3d
 nx=4, ny=4, nz=4, maximum_steps=100,
 x_upper=1e-4, y_upper=1e-4, z_upper=1e-4,
 plane_position=3.95e-5, final_time=2e-9, cfl=.01,
 riemann_solver='pelec', hydro_method='plm_ctu',
 regular_density=.31, cut_density=.93,
 initial_velocity_y=2., initial_velocity_z=-1.,
 redistribution='state_redist',
"""

    if args.inflow:
        base += """boundary_faces='outflow','outflow','inflow','outflow','outflow','outflow',
 inflow_density(3)=.31, inflow_pressure(3)=135000., inflow_velocity(2,3)=200.,
"""

    if args.selected:
        base += """thermo_model='selected', composition_count=7,
 composition_species='N2','H2O','OH','O2','O','H','H2',
 composition_mole_fractions=.55643,0.,1e-5,.14784,1e-5,1e-5,.29570,
"""

    def run(name, extra='', success=True):
        path = work / f'{name}.nml'
        path.write_text(base + f"output_file='{name}.csv',\n" + extra + '\n/\n')
        result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                capture_output=True, text=True, timeout=90)
        log = result.stdout + result.stderr
        (work / f'{name}.log').write_text(log)
        assert (result.returncode == 0) == success, (name, log)
        if not success:
            assert not (work / f'{name}.csv').exists()
            return log
        with (work / f'{name}.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 64
        for row in rows:
            values = [float(v) for v in row.values()]
            assert all(math.isfinite(v) for v in values)
            assert abs(sum(float(v) for k, v in row.items() if k.startswith('Y_')) - 1) < 5e-12
        ledgers = [line.split(':', 1)[1].split() for line in log.splitlines()
                   if line.startswith('EB inventory initial/final/boundary/wall/chemical/residual:')]
        assert ledgers
        for ledger in ledgers:
            values = [float(v) for v in ledger[1:]]
            scale = max(1e-30, *(abs(v) for v in values[:5]))
            assert abs(values[-1]) < 2e-10 * scale + 1e-25, (name, ledger)
        return (work / f'{name}.csv').read_bytes()

    inert = run('inert')
    reacting = run('reacting', 'chemistry_enabled=.true.')
    assert reacting != inert
    coupled = 'chemistry_enabled=.true., transport_enabled=.true., '
    reference = run('reference', coupled)
    assert reference != reacting
    stopped = run('stopped', coupled + "checkpoint_interval_steps=1, stop_after_checkpoint=.true., checkpoint_file='step.chk'")
    assert stopped != reference
    assert run('resumed', coupled + "restart_file='step.chk'") == reference
    assert ('Invalid EB boundary configuration' if args.inflow else 'Failed to restart public 3D EB state') in run(
        'method_mismatch', coupled + "hydro_method='pcm', restart_file='step.chk'", False)
    assert ('Invalid EB boundary configuration' if args.inflow else 'requires the pelec Riemann solver') in run('solver_mismatch', "riemann_solver='hllc'", False)
    assert ('Invalid EB boundary configuration' if args.inflow else 'hydro_method must be') in run('invalid_method', "hydro_method='unknown'", False)
    if args.inflow:
        for name, changed in [('inlet_mismatch', 'inflow_velocity(2,3)=201.'),
                              ('boundary_mismatch', "boundary_faces(3)='outflow'")]:
            assert 'Failed to restart public 3D EB state' in run(
                name, coupled + changed + ", restart_file='step.chk'", False)
        assert 'Invalid EB boundary configuration' in run('invalid_inlet', 'inflow_density(3)=-1.', False)
    print(f'PASS: EB CTU hydro, reacting/coupled stepping, exact restart and input rejection; artifacts={work}')


if __name__ == '__main__':
    main()
