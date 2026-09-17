#!/usr/bin/env python3
"""Check temperature-triggered EB layouts, composite budgets and restart."""
import argparse
import csv
import math
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser()
    for name in ('exe', 'case', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    exe, work = args.exe.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    source = args.case.read_text()
    output, budget = work / 'eb_adaptive.csv', work / 'eb_adaptive.csv.budget.csv'

    def options(text, extra):
        return text.replace('\n/', '\n ' + extra + '\n/')

    def run(name, text, success=True):
        (work / name).write_text(text)
        result = subprocess.run([str(exe), name], cwd=work, capture_output=True, text=True, timeout=90)
        log = result.stdout + result.stderr
        (work / (name + '.log')).write_text(log)
        assert (result.returncode == 0) == success, log
        return log

    for subcycle in (False, True):
        for hot in (True, False):
            name = f'{int(subcycle)}{int(hot)}'
            case = source if hot else source.replace('tag_temperature=1000.0', 'tag_temperature=1.0e9')
            case = options(case, f'automatic_timestep={str(subcycle).lower()}, subcycling={str(subcycle).lower()},')
            log = run(name + 'full.nml', case)
            counts = [int(value) for value in re.findall(r'Temperature regrid event \d+: patches=(\d+)', log)]
            assert len(counts) == 2 and all((count > 0) == hot for count in counts)
            saved = output.read_bytes(), budget.read_bytes()
            rows = list(csv.DictReader(output.read_text().splitlines()))
            assert any(int(row['level']) == 1 for row in rows) == hot
            if not hot:
                assert len(rows) == 48
            fields = [key for key in rows[0] if key.startswith('rho')]
            totals = [0.0] * len(fields)
            for row in rows:
                assert all(math.isfinite(float(value)) for value in row.values())
                assert float(row['rho']) > 0 and float(row['temperature']) > 0
                volume = float(row['kappa']) / (4 if int(row['level']) == 0 else 8) ** 3
                for i, key in enumerate(fields):
                    totals[i] += volume * float(row[key])
            entries = list(csv.DictReader(budget.read_text().splitlines()))
            assert len(entries) == len(fields)
            for total, row in zip(totals, entries, strict=True):
                values = {key: float(row[key]) for key in ('initial', 'final', 'boundary', 'wall', 'chemistry')}
                assert math.isclose(total, values['final'], rel_tol=3e-12, abs_tol=3e-12)
                residual = values['final'] - values['initial'] - sum(values[key] for key in ('boundary', 'wall', 'chemistry'))
                assert abs(residual) <= 5e-10 * max(1., *(abs(value) for value in values.values()))
            assert any(abs(float(row['chemistry'])) > 1e-15 for row in entries[5:])
            for interval in (1, 2):
                run(name + str(interval) + 'stop.nml', options(case, f"checkpoint_file='adaptive.chk', checkpoint_interval_steps={interval}, stop_after_checkpoint=.true.,"))
                checkpoint = work / 'adaptive.chk'
                original = checkpoint.read_text()
                assert original.startswith('PELEF_EB_TEMPERATURE_CHECKPOINT_')
                resume = options(case, "restart_file='adaptive.chk',")
                run(name + str(interval) + 'resume.nml', resume)
                assert saved == (output.read_bytes(), budget.read_bytes()), 'adaptive restart mismatch'
                if interval == 2:
                    wrong = resume.replace('tag_buffer=0', 'tag_buffer=1')
                    run(name + 'wrong_policy.nml', wrong, False)
                    assert saved == (output.read_bytes(), budget.read_bytes())
                    assert checkpoint.read_text() == original
            if hot:
                limited = case.replace('tag_maximum_patches=16', 'tag_maximum_patches=1')
                run(name + 'capacity.nml', limited, False)
                assert saved == (output.read_bytes(), budget.read_bytes())
    print('PASS: temperature-selected EB topology, physical inventories and exact adaptive restart')


if __name__ == '__main__':
    main()
