#!/usr/bin/env python3
"""Check serial AMR physical-boundary inventories and exact restart output."""
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
        result = subprocess.run([str(executable), str(path)], cwd=work,
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
        assert ledger(reference) == ledger(resumed)
        assert any(float(value) != 0 for value in ledger(reference).split(':')[1].split())
        changed = (source.replace('1300', '1400') if label == 'inlet' else
                   source.replace('600, 600', '650, 600', 1))
        run(f'{label}_changed', changed, f"restart_file='{label}.chk',", success=False)
    print(f'artifacts={work}')
    print('amr_physical_boundary_restart=PASS')


if __name__ == '__main__':
    main()
