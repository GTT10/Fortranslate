"""Regular MOL input, reacting iteration order, wall/LES ledgers and rejection."""
import argparse
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'source', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='mol_regular_', dir=args.work.resolve()))
    text = (args.source / 'cases/reactive_boundaries_3d/sdc_reacting.nml').read_text()
    text = text.replace("time_integrator='sdc'", "time_integrator='mol'")
    text = text.replace("reconstruction='characteristic_ppm'", "reconstruction='characteristic_plm'")
    text = text.replace('final_time=1e-8', 'final_time=1e-10')
    env = dict(os.environ)
    for name in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY', 'XAUTHORITY'):
        env.pop(name, None)
    maximum = 0.0

    def run(label, content, error=None):
        nonlocal maximum
        directory = work / label
        directory.mkdir()
        inp = directory / 'input.nml'
        inp.write_text(content)
        result = subprocess.run([str(args.exe.resolve()), str(inp)], cwd=directory, env=env,
                                capture_output=True, text=True, timeout=120)
        log = result.stdout + result.stderr
        (directory / 'run.log').write_text(log)
        output = directory / 'sdc_reacting.csv'
        if error:
            assert result.returncode != 0 and error in log, log
            assert not output.exists()
            return []
        assert result.returncode == 0, log
        assert 'Time integrator: mol' in log
        residual = float(re.search(r'Maximum normalized Euler conserved change:\s+(\S+)', log).group(1))
        assert residual < 3e-12, residual
        maximum = max(maximum, residual)
        with output.open() as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 64
        for row in rows:
            assert all(math.isfinite(float(value)) for value in row.values())
        return rows

    for iterations in (1, 3):
        result = run(f'reacting_{iterations}', text.replace('sdc_iterations=2', f'sdc_iterations={iterations}'))
        assert max(abs(float(row['heat_release_W_m3'])) for row in result) > 0
    inert = text.replace('chemistry_enabled=.true.', 'chemistry_enabled=.false.')
    walls = inert.replace("boundary_condition='periodic'", "boundary_condition='periodic',"
                          "boundary_faces='no_slip_wall','no_slip_wall',"
                          "'isothermal_wall','isothermal_wall','periodic','periodic',wall_temperature=6*900")
    les = "\n&reactive_les enabled=.true.,model='dynamic',test_model='gaussian',test_ratio=4,grid_ratio=2 /\n"
    result = run('wall_les_lf', walls.replace('sdc_iterations=2', 'sdc_iterations=3,mol_lax_friedrichs=.true.') + les)
    assert all(float(row['heat_release_W_m3']) == 0 for row in result)
    # Exercise later reacting passes with a nonzero LES wall flux. They omit
    # the first new auxiliary family; the reported boundary ledger must agree.
    result = run('reacting_wall_les', walls.replace('chemistry_enabled=.false.', 'chemistry_enabled=.true.') + les)
    assert max(abs(float(row['heat_release_W_m3'])) for row in result) > 0
    run('bad_reconstruction', text.replace("reconstruction='characteristic_plm'", "reconstruction='characteristic_ppm'"),
        'MOL requires characteristic PLM/PeleC')
    run('bad_lf_mode', text.replace("time_integrator='mol'", "time_integrator='sdc',mol_lax_friedrichs=.true."),
        'MOL LF requires MOL integration')
    print(f'PASS regular MOL input/reaction/LES ledgers; maximum residual={maximum:.9e}; {work}')


if __name__ == '__main__':
    main()
