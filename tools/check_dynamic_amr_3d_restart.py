#!/usr/bin/env python3
"""Exercise a resized-patch restart and reject changed tagging context."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from compare_amr_reactive_3d_restart import require_equal


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--mpiexec', type=Path)
    parser.add_argument('--numproc-flag', default='-n')
    parser.add_argument('--mpiexec-preflag', action='append', default=[])
    parser.add_argument('--mpiexec-postflag', action='append', default=[])
    parser.add_argument('--ranks', type=int, default=4)
    parser.add_argument('--restart-ranks', type=int, default=3)
    args = parser.parse_args()
    executable = args.executable.resolve()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    template = (Path(__file__).resolve().parents[1] /
                'cases/amr_reactive_3d/dynamic_hotspot.nml').read_text()

    def run(name: str, extra: str = '', threshold: int = 1080,
            success: bool = True) -> str:
        config = template.replace("'dynamic_coarse.csv'", f"'{name}_coarse.csv'")
        config = config.replace("'dynamic_fine.csv'", f"'{name}_fine.csv'")
        config = config.replace('threshold=1080', f'threshold={threshold}')
        config = config.replace('&amr_reactive_3d', '&amr_reactive_3d\n' + extra)
        path = work / f'{name}.nml'
        path.write_text(config)
        command = [str(executable), str(path)]
        if args.mpiexec:
            ranks = 1 if name == 'reference' else args.ranks if name == 'stopped' else args.restart_ranks
            command = ([str(args.mpiexec.resolve()), args.numproc_flag, str(ranks)]
                       + args.mpiexec_preflag + command + args.mpiexec_postflag)
        result = subprocess.run(command, cwd=work,
                                capture_output=True, text=True, timeout=240)
        output = result.stdout + result.stderr
        (work / f'{name}.log').write_text(output)
        if (result.returncode == 0) != success:
            raise RuntimeError(f'{name}: unexpected exit {result.returncode}\n{output}')
        return output

    reference = run('reference')
    stopped = run('stopped', "checkpoint_interval_steps=5, checkpoint_file='dynamic.chk', "
                  'stop_after_checkpoint=.true.,')
    assert 'Regrid at step 4:' in stopped, stopped
    resumed = run('resumed', "restart_file='dynamic.chk',")
    for level in ('coarse', 'fine'):
        require_equal(work / f'reference_{level}.csv', work / f'resumed_{level}.csv', level)
    for output in (reference, resumed):
        assert 'Completed regrids: 2' in output, output
    bad = run('changed_context', "restart_file='dynamic.chk',", threshold=1081, success=False)
    assert 'Invalid dynamic 3D AMR checkpoint header' in bad, bad
    checkpoint = (work / 'dynamic.chk').read_text()
    (work / 'truncated.chk').write_text(checkpoint[:len(checkpoint) // 2])
    run('truncated', "restart_file='truncated.chk',", success=False)
    print('PASS: resized-patch exact restart, restored regrid count, context/truncation rejection')


if __name__ == '__main__':
    main()
