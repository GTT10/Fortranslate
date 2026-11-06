"""Check automatic root intervals and serial to MPI continuation after regrid."""
import argparse
import csv
import sys
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('serial', 'mpi', 'mpiexec', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='automatic_timestep_', dir=args.work.resolve()))
    results = []
    base = """&eb_tree
 root_shape=4,4,4,child_count=1,parent=1,ratio=2,
 lower(:,1)=1,1,1,upper(:,1)=2,4,4,
 geometry_type='plane',plane_axis='z',plane_position=-10,
 time_integrator='sdc',mol_mode=.true.,sdc_iterations=3,redistribution='none',
 hydro_method='plm_ctu',riemann_solver='pelec',boundary_faces=6*'periodic',
 chemistry_enabled=.false.,transport_enabled=.false.,
 automatic_timestep=.true.,initial_timestep_factor=.5,maximum_timestep_growth=1.1,
 time_step=1e-9,final_time=FINAL,maximum_steps=10,output_file='result.csv',
 CONTROL
/
"""

    def run(label, content, parallel=False, error=None):
        folder = work / label
        folder.mkdir()
        (folder / 'input.nml').write_text(content)
        command = ([str(args.mpiexec.resolve()), '-n', '2', str(args.mpi.resolve())]
                   if parallel else [str(args.serial.resolve())]) + ['input.nml']
        result = subprocess.run(command, cwd=folder, capture_output=True, text=True, timeout=120)
        log = result.stdout + result.stderr
        (folder / 'run.log').write_text(log)
        results.append({'label': label, 'command': command, 'returncode': result.returncode})
        if error:
            assert result.returncode != 0 and error in log, log
            assert not (folder / 'result.csv').exists(), 'Rejected restart published output'
        else:
            assert result.returncode == 0, log
        return folder, log

    def history(path):
        lines = path.read_text().splitlines()
        offset = lines.index('PELEF_EB_TREE_AUTOMATIC_TIMESTEP_V1')
        fields = lines[offset + 1].split()
        return [float(value) for value in fields[:3]], fields[3] == 'T'

    for regrid in (False, True):
        name = 'regrid' if regrid else 'static'
        final = '1.55e-9' if regrid else '1.655e-9'
        scheduled = ("regrid_count=1,regrid_times(1)=5e-10,regrid_child_count(1)=1,"
                     "regrid_parent(1,1)=1,regrid_ratio(1,1)=2,"
                     "regrid_lower(:,1,1)=1,1,1,regrid_upper(:,1,1)=2,4,4," if regrid else '')
        common = base.replace('FINAL', final)
        a, log = run(name + '_full', common.replace('CONTROL', scheduled +
                    "checkpoint_file='full.chk',checkpoint_interval_steps=3,"))
        assert 'Completed steps: 3 ' in log, log
        controls, post = history(a / 'full.chk')
        expected = 5.5e-10 if regrid else 6.05e-10
        assert math.isclose(controls[2], expected, rel_tol=3e-15), controls
        assert not post
        stopped, log = run(name + '_stop', common.replace('CONTROL', scheduled +
                           "checkpoint_file='stop.chk',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,"))
        controls, post = history(stopped / 'stop.chk')
        assert controls == [.5, 1.1, 5e-10] and post == regrid, (controls, post)
        restart = f"restart_file='{stopped / 'stop.chk'}',"
        b, log = run(name + '_resume', common.replace('CONTROL', scheduled + restart), parallel=True)
        assert 'Completed steps: 3 ' in log, log
        assert (a / 'result.csv').read_bytes() == (b / 'result.csv').read_bytes()
        left = list(csv.DictReader((a / 'result.csv.budget.csv').read_text().splitlines()))
        right = list(csv.DictReader((b / 'result.csv.budget.csv').read_text().splitlines()))
        assert len(left) == len(right)
        for x, y in zip(left, right):
            # MPI sums group unchanged cells differently. State rows are exact;
            # inventory reduction permits only its machine roundoff envelope.
            for key in ('component', 'initial', 'boundary', 'wall', 'chemistry', 'forcing'):
                assert x[key] == y[key], (key, x, y)
            scale = max(abs(float(x['initial'])), abs(float(x['final'])), 1e-300)
            for key in ('final', 'residual'):
                assert abs(float(x[key]) - float(y[key])) <= 128 * sys.float_info.epsilon * scale, (key, x, y)
        run(name + '_changed_growth', common.replace('maximum_timestep_growth=1.1',
             'maximum_timestep_growth=1.2').replace('CONTROL', scheduled + restart), parallel=True,
             error='Restart automatic timestep controls mismatch')
        text = (stopped / 'stop.chk').read_text()
        lines = text.splitlines(keepends=True)
        offset = lines.index('PELEF_EB_TREE_AUTOMATIC_TIMESTEP_V1\n')
        old = work / (name + '_old.chk')
        old.write_text(''.join(lines[:offset] + lines[offset + 2:]))
        run(name + '_missing_history', common.replace('CONTROL', scheduled + f"restart_file='{old}',"),
            error='Restart automatic timestep history missing')
    # Native stop-time epsilon can finish slightly beyond the usual max_dt cap.
    content = base.replace('FINAL', '1.0001e-9').replace('initial_timestep_factor=.5',
                           'initial_timestep_factor=1').replace('CONTROL',
                           "checkpoint_file='epsilon.chk',checkpoint_interval_steps=1,")
    epsilon, log = run('stop_epsilon', content)
    assert 'Completed steps: 1 ' in log, log
    controls, post = history(epsilon / 'epsilon.chk')
    assert controls[2] == 1.0001e-9 and not post
    resumed, log = run('stop_epsilon_resume', content.replace(
        "checkpoint_file='epsilon.chk',checkpoint_interval_steps=1,",
        f"restart_file='{epsilon / 'epsilon.chk'}',"), parallel=True)
    assert (epsilon / 'result.csv').read_bytes() == (resumed / 'result.csv').read_bytes()
    (work / 'checks.json').write_text(json.dumps(results, indent=2) + '\n')
    print('PASS: automatic timestep growth, regrid hold and 1-to-2-rank continuation', work)


if __name__ == '__main__':
    main()
