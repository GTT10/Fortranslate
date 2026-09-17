#!/usr/bin/env python3
"""Check event-aligned EB remapping and restart before/at/after events."""
import argparse
import csv
import math
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    for key in ('exe', 'case', 'work'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    exe = args.exe.resolve()
    source = args.case.read_text()
    output = work / 'eb_scheduled.csv'
    budget = work / 'eb_scheduled.csv.budget.csv'

    def run(name, text, success=True):
        (work / name).write_text(text)
        result = subprocess.run([str(exe), name], cwd=work, capture_output=True, text=True, timeout=90)
        log = result.stdout + result.stderr
        (work / (name + '.log')).write_text(log)
        assert (result.returncode == 0) == success, log
        return log

    def options(text, extra):
        return text.replace('\n/', '\n ' + extra + '\n/')

    for mode in (False, True):
        case = options(source, f'automatic_timestep={str(mode).lower()}, subcycling={str(mode).lower()},')
        prefix = str(int(mode))
        log = run(prefix + 'full.nml', case)
        assert 'Completed regrids: 2' in log
        saved = output.read_bytes(), budget.read_bytes()
        rows = list(csv.DictReader(output.read_text().splitlines()))
        assert len(rows) == 160, 'final scheduled topology cell count'
        for patch in (1, 2):
            fine = [row for row in rows if int(row['patch']) == patch]
            assert len(fine) == 64, 'final fine patch size'
            assert all(int(row['level']) == 1 for row in fine)
        entries = list(csv.DictReader(budget.read_text().splitlines()))
        for row in entries:
            values = {key: float(row[key]) for key in ('initial', 'final', 'boundary', 'wall', 'chemistry')}
            assert all(math.isfinite(value) for value in values.values())
            residual = values['final'] - values['initial'] - sum(values[key] for key in ('boundary', 'wall', 'chemistry'))
            assert abs(residual) <= 5e-10 * max(1., *(abs(value) for value in values.values()))
        assert any(abs(float(row['chemistry'])) > 1e-15 for row in entries[5:])
        for interval in (1, 2, 3):
            name = prefix + str(interval)
            run(name + 'stop.nml', options(case, f"checkpoint_file='schedule.chk', checkpoint_interval_steps={interval}, stop_after_checkpoint=.true.,"))
            checkpoint = work / 'schedule.chk'
            original = checkpoint.read_text()
            lines = original.splitlines()
            assert lines[1] == 'REGRID_SCHEDULE'
            assert int(lines[2].split()[1]) == interval - 1
            resume = options(case, "restart_file='schedule.chk',")
            run(name + 'resume.nml', resume)
            assert saved == (output.read_bytes(), budget.read_bytes()), 'restart state/budget mismatch'
            if interval == 2:
                changed = resume.replace('1.5e-8', '1.6e-8')
                run(name + 'wrong_schedule.nml', changed, False)
                assert saved == (output.read_bytes(), budget.read_bytes())
                header = lines[2].split()
                header[1] = '0'
                lines[2] = ' '.join(header)
                checkpoint.write_text('\n'.join(lines) + '\n')
                run(name + 'wrong_stage.nml', resume, False)
                assert saved == (output.read_bytes(), budget.read_bytes())
                checkpoint.write_text(original)
    print('PASS: scheduled EB remap, budget, event restart and policy rejection')


if __name__ == '__main__':
    main()
