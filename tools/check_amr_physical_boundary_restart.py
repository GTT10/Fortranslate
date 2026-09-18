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
    parser.add_argument('--reacting', action='store_true')
    parser.add_argument('--pressure-outlet', action='store_true')
    parser.add_argument('--touching', action='store_true')
    parser.add_argument('--dynamic-touching', action='store_true')
    parser.add_argument('--load-balance', action='store_true')
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

    labels = ('dynamic', 'elementary') if args.reacting else ('wall', 'dynamic', 'inlet')
    if args.touching:
        labels = ('wall', 'elementary') if args.reacting else ('wall', 'inlet')
    if args.pressure_outlet:
        labels = ('wall',)
    for label in labels:
        source = template
        if args.reacting:
            source = source.replace('chemistry_enabled = .false.', 'chemistry_enabled = .true.')
            source = source.replace('final_time = 1.0e-8', 'final_time = 1.0e-10')
            source = source.replace('cfl = 0.02', 'cfl = 0.0002')
        if label == 'elementary':
            source = source.replace('full_h2o2', 'elementary')
            source = source.replace('transport_enabled = .true.', 'transport_enabled = .false.')
            source = source.replace('wall_velocity(2,1) = 20', 'wall_velocity(2,1) = 0')
            source = source.replace('isothermal_wall', 'no_slip_wall')
        if label == 'dynamic':
            for axis in ('nx', 'ny', 'nz'):
                source = source.replace(f'{axis} = 4', f'{axis} = 6')
            source = source.replace('&amr_reactive_3d', '&amr_reactive_3d\n'
                                    'regrid_interval_steps=1, regrid_temperature_threshold=1, regrid_buffer_cells=0,')
        if label == 'inlet':
            source = source.replace("'isothermal_wall'", "'inflow'")
            source = source.replace('wall_velocity(2,1) = 20',
                                    'inflow_temperature(1)=1300, inflow_velocity(1,1)=30')
        if args.touching:
            source = source.replace('coarse_i_lower = 2', 'coarse_i_lower = 1')
            if label == 'elementary':
                source = source.replace("'periodic'", "'no_slip_wall'")
                source = source.replace('coarse_j_lower = 2', 'coarse_j_lower = 1')
                source = source.replace('coarse_k_lower = 2', 'coarse_k_lower = 1')
        if args.dynamic_touching:
            source = source.replace('&amr_reactive_3d', '&amr_reactive_3d\n'
                                    'regrid_physical_boundaries=.true., regrid_interval_steps=1, regrid_temperature_threshold=1, regrid_buffer_cells=0,')
        if args.load_balance:
            source = source.replace('nx = 4', 'nx = 8')
            source = source.replace('&amr_reactive_3d', '&amr_reactive_3d\nload_balance=.true.,')
        if args.pressure_outlet:
            source = source.replace('&reactive_3d', '&reactive_3d\n outlet_pressure(2)=102325,')
            source = source.replace('coarse_i_upper = 3', 'coarse_i_upper = 4')
            source = source.replace('cfl = 0.02', 'cfl = 0.005')
            source = source.replace('&amr_reactive_3d', '&amr_reactive_3d\n'
                                    'regrid_physical_boundaries=.true., regrid_interval_steps=1, regrid_temperature_threshold=1, regrid_buffer_cells=0,')
        reference = run(f'{label}_reference', source)
        stopped = run(f'{label}_stopped', source,
                      f"checkpoint_interval_steps=1, checkpoint_file='{label}.chk', stop_after_checkpoint=.true.,")
        assert 'Stopped after checkpoint:' in stopped
        if args.load_balance:
            assert 'Balanced coarse columns:' in stopped
            if args.dynamic_touching:
                assert stopped.count('Balanced coarse columns:') >= 2
        if label == 'dynamic' or args.dynamic_touching or args.pressure_outlet:
            assert 'Regrid at step 0:' in stopped
        resumed = run(f'{label}_resumed', source, f"restart_file='{label}.chk',")
        if args.pressure_outlet:
            assert 'PRESSURE_OUTLET_V1' in (work / f'{label}.chk').read_text()
            steps = lambda log: int(next(line.split(':')[1] for line in log.splitlines()
                                        if line.startswith('Completed coarse steps:')))
            assert steps(reference) > 1
            for suffix, changed_pressure in [('pressure_changed', '103325'), ('pressure_removed', '0')]:
                changed_source = source.replace('outlet_pressure(2)=102325',
                                                f'outlet_pressure(2)={changed_pressure}')
                run(f'{label}_{suffix}', changed_source, f"restart_file='{label}.chk',", success=False)

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
        if args.dynamic_touching:
            checkpoint = work / f'{label}.chk'
            original_checkpoint = checkpoint.read_text()
            assert ' PHYSICAL_TAGGING_V1' in original_checkpoint
            try:
                checkpoint.write_text(original_checkpoint.replace(' PHYSICAL_TAGGING_V1', '', 1))
                run(f'{label}_tagging_policy', source, f"restart_file='{label}.chk',", success=False)
            finally:
                checkpoint.write_text(original_checkpoint)
        changed = (source.replace('1300', '1400') if label == 'inlet' else
                   source.replace('600, 600', '650, 600', 1))
        run(f'{label}_changed', changed, f"restart_file='{label}.chk',", success=False)
        if args.reacting:
            inactive = source.replace('chemistry_enabled = .true.', 'chemistry_enabled = .false.')
            run(f'{label}_inactive_reference', inactive)
            assert (work / f'{label}_inactive_reference_coarse.csv').read_bytes() != (work / f'{label}_reference_coarse.csv').read_bytes()
            changed_tolerance = source.replace('2.0e-7', '3.0e-7')
            run(f'{label}_tolerance', changed_tolerance, f"restart_file='{label}.chk',", success=False)
            checkpoint = work / f'{label}.chk'
            original_checkpoint = checkpoint.read_text()
            assert 'FIXED_CONTEXT' in original_checkpoint
            wrong_integrator = 'implicit' if label == 'elementary' else 'explicit'
            lines = original_checkpoint.splitlines()
            context_line = lines.index('FIXED_CONTEXT')
            lines[context_line + 2] = wrong_integrator
            checkpoint.write_text('\n'.join(lines) + '\n')
            run(f'{label}_integrator', source, f"restart_file='{label}.chk',", success=False)
            lines = original_checkpoint.splitlines()
            reaction_rates = lines.index('REACTIONS') + 6
            rates = lines[reaction_rates].split()
            rates[0] = str(float(rates[0]) * 1.01)
            lines[reaction_rates] = ' '.join(rates)
            checkpoint.write_text('\n'.join(lines) + '\n')
            run(f'{label}_reaction', source, f"restart_file='{label}.chk',", success=False)
            lines = original_checkpoint.splitlines()
            identity = lines[context_line + 1]
            lines[context_line + 1] = identity[:-1] + ('0' if identity[-1] != '0' else '1')
            checkpoint.write_text('\n'.join(lines) + '\n')
            run(f'{label}_identity', source, f"restart_file='{label}.chk',", success=False)
            checkpoint.write_text(original_checkpoint)

    print(f'artifacts={work}')
    print('amr_physical_boundary_restart=PASS')


if __name__ == '__main__':
    main()
