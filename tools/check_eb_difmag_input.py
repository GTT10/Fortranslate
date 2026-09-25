#!/usr/bin/env python3
"""Check EB coefficient input and exact checkpoint continuation."""
import argparse
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--tree', action='store_true')
    parser.add_argument('--two-level', choices=('reference', 'multipatch'))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = (root / 'cases/reactive_eb_3d/density_sheet_state.nml').read_text()
    source = source.replace('maximum_steps = 10', 'maximum_steps = 100')
    source = source.replace('cfl = 0.5', 'cfl = 0.01').replace('final_time = 5.0e-5', 'final_time = 1.0e-5')
    if args.tree:
        source = (root / 'cases/reactive_eb_tree_3d/input.nml').read_text()
        source = source.replace('&eb_tree', '&eb_tree\n subcycling=.true.,', 1)
        source = source.replace("output_file='eb_tree.csv'", "output_file = 'reactive_eb_3d_state.csv'")
    if args.two_level:
        source = (root / f'cases/reactive_eb_two_level_3d/{args.two_level}.nml').read_text()
        source = source.replace('&eb_two_level', '&eb_two_level\n subcycling=.true.,', 1)
        source = source.replace('automatic_timestep=.true.', 'automatic_timestep=.false.')
        for output in ('eb_two_level.csv', 'eb_multipatch.csv'):
            source = source.replace(f"output_file='{output}'", "output_file = 'reactive_eb_3d_state.csv'")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='run-', dir=args.work_dir.resolve()))

    def run(name, extra='', success=True):
        text = source.replace("output_file = 'reactive_eb_3d_state.csv'", f"output_file = '{name}.csv'")
        text = text.rsplit('/', 1)[0] + extra + '\n/\n'
        path = work / f'{name}.nml'
        path.write_text(text)
        result = subprocess.run([str(args.executable.resolve()), str(path)], cwd=work,
                                text=True, capture_output=True, timeout=90)
        log = result.stdout + result.stderr
        (work / f'{name}.log').write_text(log)
        assert (result.returncode == 0) == success, (name, log)
        if not success:
            assert not (work / f'{name}.csv').exists()
            return log
        return (work / f'{name}.csv').read_bytes()

    reference = run('reference', 'difmag=.037')
    stopped = run('stopped', "difmag=.037, checkpoint_interval_steps=1, stop_after_checkpoint=.true., checkpoint_file='step.chk'")
    assert stopped != reference, 'checkpoint must stop before final output'
    assert run('resumed', "difmag=.037, restart_file='step.chk'") == reference
    log = run('changed', "difmag=.05, restart_file='step.chk'", False)
    diagnostic = 'Invalid tree checkpoint record' if args.tree else 'Failed to restart public 3D EB state'
    if args.two_level:
        diagnostic = 'Invalid two-level patch record'
    assert diagnostic in log
    for label, value in [('negative', '-1'), ('nonfinite', 'NaN')]:
        assert 'Invalid artificial-viscosity coefficient' in run(label, f'difmag={value}', False)
    baseline = run('baseline')
    assert run('zero', 'difmag=0') == baseline
    run('legacy_stop', "checkpoint_interval_steps=1, stop_after_checkpoint=.true., checkpoint_file='zero.chk'")
    assert run('zero_resume', "difmag=0, restart_file='zero.chk'") == baseline
    print(f'PASS: EB coefficient input and exact restart; artifacts={work}')


if __name__ == '__main__':
    main()
