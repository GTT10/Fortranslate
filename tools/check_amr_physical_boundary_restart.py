#!/usr/bin/env python3
"""Check serial/MPI AMR physical-boundary inventories and exact restart output."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

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
    parser.add_argument('--restart-ranks', type=int, default=2)
    parser.add_argument('--reference-executable', type=Path)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='run-', dir=args.work_dir.resolve()))
    executable = args.executable.resolve()
    template = (Path(__file__).resolve().parents[1] /
                'cases/amr_reactive_3d/physical_wall.nml').read_text()

    def run(name: str, source: str, extra: str = '', success: bool = True) -> str:
        source = source.replace('amr_physical_wall_coarse.csv', f'{name}_coarse.csv')
        source = source.replace('amr_physical_wall_fine.csv', f'{name}_fine.csv')
        source = source.replace('&amr_reactive_3d', '&amr_reactive_3d\n' + extra)
        path = work / f'{name}.nml'
        path.write_text(source)
        command = [str(executable), str(path)]
        if name.endswith('_reference') and args.reference_executable:
            command[0] = str(args.reference_executable.resolve())
        elif args.mpiexec:
            ranks = args.ranks if name.endswith('_stopped') else args.restart_ranks
            command = ([str(args.mpiexec.resolve()), args.numproc_flag, str(ranks)]
                       + args.mpiexec_preflag + command + args.mpiexec_postflag)
        result = subprocess.run(command, cwd=work,
                                capture_output=True, text=True, timeout=180)
        output = result.stdout + result.stderr
        (work / f'{name}.log').write_text(output)
        if (result.returncode == 0) != success:
            raise AssertionError(f'{name}: exit {result.returncode}\n{output}')
        if not success:
            assert not (work / f'{name}_coarse.csv').exists()
        return output

    for label in ('wall', 'dynamic', 'inlet'):
        source = template
        if label == 'dynamic':
            for axis in ('nx', 'ny', 'nz'):
                source = source.replace(f'{axis} = 4', f'{axis} = 6')
            source = source.replace('&amr_reactive_3d', '&amr_reactive_3d\n'
                                    'regrid_interval_steps=1, regrid_temperature_threshold=1, regrid_buffer_cells=0,')
        if label == 'inlet':
            source = source.replace("'isothermal_wall'", "'inflow'")
            source = source.replace('wall_velocity(2,1) = 20',
                                    'inflow_temperature(1)=1300, inflow_velocity(1,1)=30')
        reference = run(f'{label}_reference', source)
        stopped = run(f'{label}_stopped', source,
                      f"checkpoint_interval_steps=1, checkpoint_file='{label}.chk', stop_after_checkpoint=.true.,")
        assert 'Stopped after checkpoint:' in stopped
        if label == 'dynamic':
            assert 'Regrid at step 0:' in stopped
        resumed = run(f'{label}_resumed', source, f"restart_file='{label}.chk',")
        for level in ('coarse', 'fine'):
            require_equal(work / f'{label}_reference_{level}.csv',
                          work / f'{label}_resumed_{level}.csv', f'{label}_{level}')
        ledger = lambda log: next(line for line in log.splitlines()
                                  if line.startswith('Conserved boundary change:'))
        if args.mpiexec:
            original = [float(value) for value in ledger(reference).split(':')[1].split()]
            restored = [float(value) for value in ledger(resumed).split(':')[1].split()]
            assert len(original) == len(restored)
            scales = [max(abs(a), abs(b)) for a, b in zip(original, restored)]
            # Momentum is a vector; use a shared scale for cancelling components.
            scales[1:4] = [max(scales[1:4])] * 3
            assert all(abs(a-b) <= 2e-11 * scale + 1e-25
                       for a, b, scale in zip(original, restored, scales)), (original, restored)
        else:
            assert ledger(reference) == ledger(resumed)
        assert any(float(value) != 0 for value in ledger(reference).split(':')[1].split())
        changed = (source.replace('1300', '1400') if label == 'inlet' else
                   source.replace('600, 600', '650, 600', 1))
        run(f'{label}_changed', changed, f"restart_file='{label}.chk',", success=False)
    print(f'artifacts={work}')
    print('amr_physical_boundary_restart=PASS')


if __name__ == '__main__':
    main()
