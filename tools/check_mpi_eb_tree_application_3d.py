#!/usr/bin/env python3
"""Serial/MPI EB tree state parity, rank-changing restart and output protection."""
import argparse
import csv
import math
import os
import re
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial', type=Path, required=True)
    parser.add_argument('--mpi', type=Path, required=True)
    parser.add_argument('--mpiexec', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--automatic', action='store_true')
    parser.add_argument('--scheduled', action='store_true')
    parser.add_argument('--adaptive', action='store_true')
    parser.add_argument('--subcycling', action='store_true')
    parser.add_argument('--mach', action='store_true')
    parser.add_argument('--selected', action='store_true')
    parser.add_argument('--rebalance', action='store_true')
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    template = args.input.read_text().replace("output_file='eb_tree.csv'", "output_file='result.csv'")
    if args.subcycling:
        template = template.replace('&eb_tree', '&eb_tree\n subcycling=.true.,', 1)
    if args.automatic:
        template = template.replace('&eb_tree', '&eb_tree\n automatic_timestep=.true.,hydro_cfl=5.e-5,', 1)
    if args.subcycling and args.automatic:
        template = template.replace('hydro_cfl=5.e-5', 'hydro_cfl=1.e-5')
    if args.mach:
        template = template.replace('&eb_tree', '&eb_tree\n inflow_mach(1)=.2,inflow_mach(4)=.1,'
                                    'inflow_temperature(1)=400,inflow_temperature(4)=400,'
                                    'inflow_pressure(4)=120000,', 1)
    env = os.environ.copy()
    env['OMPI_MCA_rmaps_base_oversubscribe'] = '1'

    def run(name, ranks=0, extra='', text=None, success=True):
        directory = work / name
        directory.mkdir(exist_ok=True)
        source = template if text is None else text
        if ranks and args.rebalance:
            extra += 'rebalance_interval_steps=1,rebalance_minimum_gain=0.,'
        (directory / 'input.nml').write_text(source.rsplit('/', 1)[0] + extra + '\n/\n')
        command = [str(args.serial.resolve()), 'input.nml']
        if ranks:
            command = [str(args.mpiexec.resolve()), '-n', str(ranks), str(args.mpi.resolve()), 'input.nml']
        result = subprocess.run(command, cwd=directory, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
        (directory / 'run.log').write_text(result.stdout)
        if (result.returncode == 0) != success:
            raise AssertionError(f'{name}: exit {result.returncode}\n{result.stdout}')
        return directory

    serial = run('serial')
    full = run('mpi2', 2)
    third = run('mpi3', 3)
    stop_policy = "checkpoint_file='saved.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true.,"
    stopped = run('stopped2', 2, stop_policy)
    if args.automatic:
        match = re.search(r'Minimum accepted interval:\s+(\S+)', (stopped / 'run.log').read_text())
        assert match and 0 < float(match[1]) < 1e-9, 'automatic interval did not shorten first step'
    resumed = run('resumed2', 2, f"restart_file='{stopped / 'saved.chk'}',")
    changed = run('resumed3', 3, f"restart_file='{stopped / 'saved.chk'}',")
    to_serial = run('mpi_to_serial', 0, f"restart_file='{stopped / 'saved.chk'}',")
    serial_stop = run('serial_stop', 0, stop_policy)
    from_serial = run('serial_to_mpi', 3, f"restart_file='{serial_stop / 'saved.chk'}',")
    patch_marker = b'PELEF_REACTIVE_EB_3D_CHECKPOINT\n'
    mpi_saved = (stopped / 'saved.chk').read_bytes().split(patch_marker, 1)
    serial_saved = (serial_stop / 'saved.chk').read_bytes().split(patch_marker, 1)
    assert len(mpi_saved) == len(serial_saved) == 2
    assert mpi_saved[1] == serial_saved[1], 'checkpoint patch records differ from serial'
    expected = (serial / 'result.csv').read_bytes()
    for directory in (full, third, resumed, changed, to_serial, from_serial):
        assert (directory / 'result.csv').read_bytes() == expected, directory
    for filename in ('result.csv', 'result.csv.budget.csv'):
        if args.rebalance and filename.endswith('.budget.csv'):
            continue  # Timings can change owners and MPI reduction grouping.
        assert (full / filename).read_bytes() == (resumed / filename).read_bytes(), filename

    scheduled_results = []
    if args.scheduled or args.adaptive:
        event_nodes = ((2, 5), (3, 5)) if args.adaptive else ((2, 1), (3, 4))
        for step, nodes in event_nodes:
            stop = run(f'stop_event_{step}', 2, stop_policy.replace('steps=1', f'steps={step}'))
            assert f'nodes: {nodes}' in (stop / 'run.log').read_text()
            resume = run(f'resume_event_{step}', 3, f"restart_file='{stop / 'saved.chk'}',")
            assert (resume / 'result.csv').read_bytes() == expected, step
            scheduled_results.append(resume)

    reference_budget = list(csv.DictReader((serial / 'result.csv.budget.csv').open()))
    worst = 0.0
    for directory in (serial, full, third, resumed, changed, to_serial, from_serial, *scheduled_results):
        with (directory / 'result.csv').open() as stream:
            cells = list(csv.DictReader(stream))
        with (directory / 'result.csv.budget.csv').open() as stream:
            budget = list(csv.DictReader(stream))
        fields = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE'] + [key for key in cells[0] if key.startswith('rhoY_')]
        assert len(fields) == len(budget) == len(reference_budget)
        for field, row, reference in zip(fields, budget, reference_budget):
            initial, final = float(row['initial']), float(row['final'])
            inventory = math.fsum(float(cell[field]) * float(cell['fluid_volume']) for cell in cells)
            assert math.isfinite(inventory) and abs(inventory-final) <= 3e-14*max(1, abs(final))
            exchange = sum(float(row[key]) for key in ('boundary', 'wall', 'chemistry'))
            residual = abs(final-initial-exchange)/max(1, abs(initial), abs(final), abs(exchange))
            assert residual <= 5e-10
            worst = max(worst, residual)
            for key in ('initial', 'final', 'boundary', 'wall', 'chemistry'):
                assert abs(float(row[key])-float(reference[key])) <= 3e-14*max(1, abs(initial)), (directory, key)

    if args.adaptive:
        cold = template.replace('1000,1.e9', '1.e9,1.e9')
        assert cold != template, 'missing adaptive fixture thresholds'
        cold_serial = run('cold_serial', text=cold)
        cold_mpi = run('cold_mpi', 3, text=cold)
        assert (cold_serial / 'result.csv').read_bytes() == (cold_mpi / 'result.csv').read_bytes()
        with (cold_mpi / 'result.csv').open() as stream:
            assert {int(row['level']) for row in csv.DictReader(stream)} == {0}
        run('capacity', 2, text=template.replace('tag_maximum_nodes=33', 'tag_maximum_nodes=1'), success=False)

    broken = work / 'broken.chk'
    broken.write_text((stopped / 'saved.chk').read_text().replace('END_EB_TREE_CHECKPOINT', 'BROKEN_END'))
    damaged = work / 'damaged_patch.chk'
    prefix, suffix = (stopped / 'saved.chk').read_text().rsplit('\nEND_CHECKPOINT', 1)
    rows = prefix.splitlines()
    values = rows[-1].split()
    assert len(values) > 6, 'missing final patch state row'
    values[0] = '-1.0'
    # Keep V2 offsets valid so rejection exercises the patch state validator.
    replacement = ' '.join(values)
    assert len(replacement) <= len(rows[-1])
    rows[-1] = replacement.ljust(len(rows[-1]))
    damaged.write_text('\n'.join(rows) + '\nEND_CHECKPOINT' + suffix)
    failures = [
        ('corrupt', f"restart_file='{broken}',"),
        ('damaged_patch', f"restart_file='{damaged}',"),
        ('wrong_policy', f"restart_file='{stopped / 'saved.chk'}',time_step=2.e-9,"),
        ('wrong_auto', f"restart_file='{stopped / 'saved.chk'}',automatic_timestep={'.false.' if args.automatic else '.true.'},"),
        ('wrong_cfl', f"restart_file='{stopped / 'saved.chk'}',hydro_cfl=1.e-6,"),
        ('wrong_subcycle', f"restart_file='{stopped / 'saved.chk'}',subcycling={'.false.' if args.subcycling else '.true.'},"),
    ]
    indexed_rows = (stopped / 'saved.chk').read_text().splitlines(keepends=True)
    tree_header = indexed_rows.index('PELEF_EB_TREE_CHECKPOINT_V2\n')
    patch_count = int(indexed_rows[tree_header + 1])
    first_length = tree_header + patch_count + 4
    first_patch = first_length + patch_count
    assert indexed_rows[first_patch] == 'PELEF_REACTIVE_EB_3D_CHECKPOINT\n'
    for name, length in [('negative_index', -1), ('oversized_index', 2**63-1),
                         ('shifted_index', int(indexed_rows[first_length]) + 1),
                         ('balanced_index', int(indexed_rows[first_length]) + 1)]:
        rows = indexed_rows.copy()
        rows[first_length] = f'{length:20d}\n'
        if name == 'balanced_index':
            assert patch_count > 1
            rows[first_length + 1] = f'{int(rows[first_length + 1]) - 1:20d}\n'
        invalid_index = work / f'{name}.chk'
        invalid_index.write_text(''.join(rows))
        failures.append((name, f"restart_file='{invalid_index}',"))
        run(name + '_serial', extra=f"restart_file='{invalid_index}',", success=False)
    if args.mach:
        for field, value in [('inflow_mach(1)', '.3'), ('inflow_temperature(1)', '401'),
                             ('inflow_pressure(4)', '120001')]:
            failures.append((f'wrong_{field[:field.index("(")]}',
                             f"restart_file='{stopped / 'saved.chk'}',{field}={value},"))
        failures.append(('invalid_mach', 'inflow_mach(1)=101,'))
        failures.append(('invalid_inlet_temperature', 'inflow_temperature(1)=0,'))
    if args.selected:
        context = (stopped / 'saved.chk').read_text()
        assert 'ELEMENT_CONTEXT' in context
        match = re.search(r'SELECTED_CONTEXT\n([0-9a-f]{64})', context)
        assert match, 'selected mechanism identity missing'
        sha = match[1]
        altered = ('0' if sha[0] != '0' else '1') + sha[1:]
        bad_identity = work / 'bad_identity.chk'
        bad_identity.write_text(context.replace(sha, altered))
        failures.extend([
            ('wrong_identity', f"restart_file='{bad_identity}',"),
            ('wrong_composition', f"restart_file='{stopped / 'saved.chk'}',composition_mole_fractions=0.30,0.15,0.55,"),
        ])
    if args.adaptive:
        failures.append(('wrong_tag', f"restart_file='{stopped / 'saved.chk'}',tag_buffer=1,"))
    for name, extra in failures:
        directory = work / name
        directory.mkdir(exist_ok=True)
        paths = [directory / filename for filename in ('result.csv', 'result.csv.budget.csv', 'protected.chk')]
        for path in paths:
            path.write_text('sentinel\n')
        run(name, 2, extra + "checkpoint_file='protected.chk',checkpoint_interval_steps=1,", success=False)
        assert all(path.read_text() == 'sentinel\n' for path in paths), name
        if name == 'wrong_subcycle':
            assert 'Restart subcycling mode mismatch' in (directory / 'run.log').read_text()
        if name == 'wrong_tag':
            assert 'Tagging policy mismatch' in (directory / 'run.log').read_text()
        if name in ('wrong_auto', 'wrong_cfl'):
            expected_error = 'Restart interval mode mismatch' if name == 'wrong_auto' else 'Restart CFL mismatch'
            assert expected_error in (directory / 'run.log').read_text(), name
    print(f'PASS: serial/2/3-rank exact state, same-rank exact budget, changed-rank and serial restart; residual={worst:.17e}')


if __name__ == '__main__':
    main()
