#!/usr/bin/env python3
"""Run a bounded native-MOL thermal-contact mesh diagnostic from saved imports.

Preserve imported rho, species, held energy and temperature. Replace momentum
with rho*(10,0,0), and set total energy to held energy plus kinetic energy.
These are new diagnostic inputs, never replacements for the pinned reference.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source, target, cells):
    lines = source.read_text().splitlines()
    if lines[0] != 'PELEF_REACTING_EB_REFERENCE_V2' or int(lines[1]) != cells:
        raise ValueError('Expected a saved V2 import at the requested resolution')
    count = 5 + 63
    if int(lines[count]) != 2 * cells**3 or len(lines) != count + 1 + 2 * cells**3:
        raise ValueError('Incomplete import')
    spacing = [float(x) for x in lines[3].split()]
    if any(abs(cells * x - .01) > 1e-15 for x in spacing):
        raise ValueError('Expected the same 1 cm domain at every resolution')
    if abs(float(lines[4].split()[0]) - .0033) > 1e-15:
        raise ValueError('Expected the same embedded plane')
    temperatures = []
    seen = set()
    output = ['PELEF_REACTING_EB_REFERENCE_V4', str(cells), '1e-12 1e-18', '1',
              '4 2 1 5e-8'] + lines[3:count + 1]
    for line in lines[count + 1:]:
        row = line.split()
        index = tuple(map(int, row[:4]))
        if index in seen or index[0] not in (1, 2) or any(not 1 <= i <= cells for i in index[1:]):
            raise ValueError('Invalid or duplicate imported index')
        seen.add(index)
        q = list(map(float, row[4:]))
        if len(q) != 26 or not all(map(math.isfinite, q)):
            raise ValueError('Invalid LiDryer state record')
        rho, held = q[0], q[14]
        if rho <= 0 or abs(sum(q[5:14]) - rho) > 1e-12 * rho:
            raise ValueError('Invalid material inventory')
        q[1:4] = [rho * 10., 0., 0.]
        q[4] = held + .5 * rho * 100.
        temperatures.append(q[15])
        output.append(' '.join(map(str, index)) + ' ' + ' '.join(f'{v:.17e}' for v in q))
    target.write_text('\n'.join(output) + '\n')
    return {'source': str(source), 'source_sha256': sha(source),
            'input': str(target), 'input_sha256': sha(target),
            'temperature_range_K': [min(temperatures), max(temperatures)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    for n in (8, 16, 32):
        parser.add_argument(f'--input{n}', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--maximum-grid32-seconds', type=float, default=3600.)
    args = parser.parse_args()
    out = args.work.resolve()
    out.mkdir(parents=True, exist_ok=False)
    probe = out / 'frozen-probe'
    shutil.copy2(args.probe, probe)
    record = {'status': 'preparing', 'qualification': False, 'pid': os.getpid(),
              'cwd': str(Path.cwd()), 'binary': str(probe), 'binary_sha256': sha(probe),
              'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
              'runner_sha256': sha(Path(__file__)), 'endpoint_s': 2e-7, 'dt_s': 5e-8,
              'chemistry': False, 'transport': False, 'velocity_m_s': [10, 0, 0],
              'clock': 'recursive', 'runs': [], 'maximum_grid32_seconds': args.maximum_grid32_seconds}

    def save():
        temp = out / 'record.tmp'
        temp.write_text(json.dumps(record, indent=2) + '\n')
        temp.replace(out / 'record.json')

    save()
    try:
        prepared = {}
        fits = None
        for n in (8, 16, 32):
            source = getattr(args, f'input{n}').resolve()
            source_fits = source.read_text().splitlines()[5:68]
            if fits is not None and source_fits != fits:
                raise ValueError('Transport fit identity differs between imported grids')
            fits = source_fits
            prepared[n] = prepare(source, out / f'grid{n}.dat', n)
        record['inputs'] = prepared
        save()
        for n in (8, 16, 32):
            if n == 32:
                # Fixed steps, 8x cells; retain 50% overhead margin in the pilot estimate.
                forecast = record['runs'][-1]['seconds'] * 8 * 1.5
                record['grid32_forecast_seconds'] = forecast
                if forecast > args.maximum_grid32_seconds:
                    record['status'] = 'cost_limit'; save(); return
            prefix = out / f'grid{n}'
            command = [str(probe), prepared[n]['input'], str(prefix), '--audit-native',
                       '--separate-momentum', '--pressure-evolve-reset',
                       '--imported-state', '--consistent-velocity']
            item = {'cells': n, 'argv': command, 'status': 'running'}
            record['runs'].append(item); record['status'] = 'running'; save()
            start = time.monotonic()
            with (out / f'grid{n}.log').open('w') as log:
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
                item['pid'] = process.pid; save()
                try:
                    item['returncode'] = process.wait(timeout=args.maximum_grid32_seconds)
                except subprocess.TimeoutExpired:
                    process.kill(); process.wait(); raise
            item['seconds'] = time.monotonic() - start
            if item['returncode'] != 0:
                raise RuntimeError(f'Grid {n} failed: see its log')
            for policy in ('velocity', 'paired'):
                with Path(str(prefix) + f'_{policy}_native.csv').open() as stream:
                    rows = list(csv.DictReader(stream))
                if len(rows) != 5 or abs(float(rows[-1]['time_s']) - 2e-7) > 1e-20:
                    raise ValueError('Incomplete trajectory')
                budget = max(float(r['max_budget_residual']) for r in rows)
                if not math.isfinite(budget) or budget > 1e-12:
                    raise ValueError('Conservation budget failed')
                initial = float(rows[0]['max_pressure_offset_Pa'])
                if not math.isfinite(initial) or initial > 1e-7:
                    raise ValueError('Initial state is not pressure equilibrated')
                item[policy] = {'endpoint_pressure_error_Pa': float(rows[-1]['max_pressure_offset_Pa']),
                                'initial_pressure_error_Pa': initial, 'maximum_budget': budget,
                                'numerical_energy_J': float(rows[-1]['numerical_energy_J'])}
            item['status'] = 'completed'; save()
            print(f'grid {n} completed in {item["seconds"]:.1f} s', flush=True)
        record['orders'] = {}
        for policy in ('velocity', 'paired'):
            errors = [r[policy]['endpoint_pressure_error_Pa'] for r in record['runs']]
            record['orders'][policy] = [math.log2(a / b) if a > 0 and b > 0 else None
                                        for a, b in zip(errors, errors[1:])]
        if sha(probe) != record['binary_sha256']:
            raise ValueError('Frozen executable changed')
        record['status'] = 'completed'
    except Exception as exc:
        record['status'] = 'failed'; record['error'] = str(exc); raise
    finally:
        save()


if __name__ == '__main__':
    main()
