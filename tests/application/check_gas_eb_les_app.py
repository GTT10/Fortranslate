"""Check gas-only EB-tree LES routing, budgets and changed-rank restart."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'mpiexec', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--recursive', action='store_true')
    parser.add_argument('--dynamic', action='store_true')
    parser.add_argument('--model', choices=('smagorinsky', 'wale', 'vreman', 'dynamic'))
    parser.add_argument('--serial-exe', type=Path)
    parser.add_argument('--selected', action='store_true')
    parser.add_argument('--nonmol', action='store_true', help='Use CTU with stage diffusion predictors')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='gas_les_', dir=args.work.resolve()))
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1')
    model = args.model or ('dynamic' if args.dynamic else 'smagorinsky')
    template = """&eb_tree
 root_shape=4,4,4,child_count=1,parent(1)=1,ratio(1)=2,
 lower(:,1)=1,1,1,upper(:,1)=2,4,4,
 domain_lower=0,0,0,domain_upper=.004,.004,.004,
 geometry_type='level_set',level_set_file='level_set.dat',
 time_step=1.e-8,final_time=2.e-8,maximum_steps=100,
 chemistry_enabled=.false.,transport_enabled=.true.,
 time_integrator='sdc',sdc_iterations=2,mol_mode=.true.,
 redistribution='state_redist',riemann_solver='pelec',hydro_method='plm_ctu',
 boundary_faces='inflow',5*'outflow',inflow_velocity(:,1)=20.,100.,0.,
 inflow_pressure(1)=135000.,
 output_file='result.csv',
 CONTROLS
/
"""
    if args.recursive:
        template = template.replace('CONTROLS', 'subcycling=.true.,CONTROLS')
    if args.nonmol:
        template = template.replace('mol_mode=.true.', 'mol_mode=.false.')
    if args.selected:
        template = template.replace('&eb_tree', "&eb_tree\n composition_count=1,"
            "composition_species(1)='N2',composition_mole_fractions(1)=1.,", 1)
    les = f"\n&eb_tree_les les_enabled=.true.,les_model='{model}',les_cs=.17 /\n"

    def run(label, controls='', settings=les, ranks=2, error=None, serial=False):
        directory = work / label
        directory.mkdir()
        (directory / 'level_set.dat').write_text('8 8 8\n' + '1.0\n' * 9**3)
        (directory / 'input.nml').write_text(template.replace('CONTROLS', controls) + settings)
        command = ([str(args.serial_exe.resolve()), 'input.nml'] if serial else
                   [str(args.mpiexec.resolve()), '-n', str(ranks), str(args.exe.resolve()), 'input.nml'])
        result = subprocess.run(command, cwd=directory,
                                env=env, capture_output=True, text=True, timeout=180)
        log = result.stdout + result.stderr
        (directory / 'run.log').write_text(log)
        if error:
            assert result.returncode != 0 and error in log, log
            assert not (directory / 'result.csv').exists()
        else:
            assert result.returncode == 0, log
        return directory

    def outputs(directory):
        return tuple((directory / name).read_bytes() for name in
                     ('result.csv', 'result.csv.budget.csv'))

    disabled = run('disabled', settings='')
    continuous = run('continuous')
    if model == 'smagorinsky':
        assert outputs(disabled)[0] != outputs(continuous)[0], 'LES must change the shear-flow gas state'
    with (continuous / 'result.csv.budget.csv').open() as stream:
        budget = list(csv.DictReader(stream))
    mass, energy = float(budget[0]['initial']), float(budget[4]['initial'])
    errors = []
    for row in budget:
        scale = max(abs(float(row[key])) for key in
                    ('initial', 'final', 'boundary', 'wall', 'chemistry', 'forcing'))
        if 2 <= int(row['component']) <= 4:
            scale = max(scale, math.sqrt(mass * energy))
        errors.append(abs(float(row['residual'])) / max(scale, 1.e-30))
    assert errors and all(math.isfinite(value) and value < 3.e-12 for value in errors)
    stopped = run('stopped', "checkpoint_file='saved',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,")
    restart = f"restart_file='{stopped / 'saved'}',"
    resumed = run('resumed', restart, ranks=3)
    assert outputs(resumed) == outputs(continuous), 'Changed-rank gas LES restart differs'
    run('disabled_restart', restart, settings='', error='Restart LES enable mismatch')
    run('coefficient_restart', restart, settings=les.replace('les_cs=.17', 'les_cs=.18'),
        error='Restart LES coefficients mismatch')
    run('model_restart', restart, settings=les.replace(model, 'smagorinsky' if model == 'dynamic' else 'dynamic'),
        error='Restart LES model mismatch')
    serial_error = None
    if args.serial_exe:
        serial_result = run('serial_continuous', serial=True)
        serial_stopped = run('serial_stopped',
            "checkpoint_file='saved',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,", serial=True)
        serial_resumed = run('serial_resumed', f"restart_file='{serial_stopped / 'saved'}',", serial=True)
        assert outputs(serial_result) == outputs(serial_resumed), 'Serial LES continuation differs'
        mpi_from_serial = run('mpi_from_serial', f"restart_file='{serial_stopped / 'saved'}',", ranks=3)
        serial_from_mpi = run('serial_from_mpi', restart, serial=True)
        def compare_gas(first, second):
            def rows(directory):
                with (directory / 'result.csv').open() as stream:
                    return list(csv.DictReader(stream))
            left, right = rows(first), rows(second)
            assert len(left) == len(right)
            maximum = 0.
            for a, b in zip(left, right):
                assert a.keys() == b.keys()
                for key in a:
                    x, y = float(a[key]), float(b[key])
                    error = abs(x - y) / max(1., abs(x), abs(y))
                    assert math.isfinite(error) and error < 3.e-12, (key, x, y)
                    maximum = max(maximum, error)
            return maximum
        serial_error = max(compare_gas(continuous, serial_result),
                           compare_gas(serial_result, mpi_from_serial),
                           compare_gas(continuous, serial_from_mpi))
        run('serial_coefficient_restart', restart,
            settings=les.replace('les_cs=.17', 'les_cs=.18'), serial=True,
            error='Restart LES coefficients mismatch')
    (work / 'report.json').write_text(json.dumps({'model': model, 'recursive': args.recursive,
        'selected': args.selected,
        'nonmol': args.nonmol,
        'maximum_budget_error': max(errors), 'changed_rank_restart': 'exact',
        'rejected_restart_controls': ['enabled', 'coefficient', 'model'],
        'serial_mpi_maximum_scaled_gas_error': serial_error}, indent=2) + '\n')
    print(work)


if __name__ == '__main__':
    main()
