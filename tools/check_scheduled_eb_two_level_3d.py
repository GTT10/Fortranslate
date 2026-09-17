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
    parser.add_argument("--coarse-gap", action="store_true")
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    exe = args.exe.resolve()
    source = args.case.read_text()
    if args.coarse_gap:
        source = source.replace("regrid_patch_count=1,2", "regrid_patch_count=0,2")
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
            if args.coarse_gap and interval == 2:
                coarse_rows = list(csv.DictReader(output.read_text().splitlines()))
                assert len(coarse_rows) == 48 and all(int(row['level']) == 0 for row in coarse_rows)
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
    if args.coarse_gap:
        # Static coarse-only initialization and its separate checkpoint schema.
        case = source.replace('patch_count=2,', 'patch_count=0,', 1).replace('regrid_count=2,', 'regrid_count=0,', 1)
        case = options(case, 'automatic_timestep=.true., subcycling=.true.,')
        run('coarse_full.nml', case)
        saved = output.read_bytes(), budget.read_bytes()
        rows = list(csv.DictReader(output.read_text().splitlines()))
        assert len(rows) == 48 and all(int(row['level']) == 0 for row in rows)
        run('coarse_stop.nml', options(case, "checkpoint_file='coarse.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true.,"))
        assert (work / 'coarse.chk').read_text().startswith('PELEF_EB_COARSE_CHECKPOINT_V4')
        run('coarse_resume.nml', options(case, "restart_file='coarse.chk',"))
        assert saved == (output.read_bytes(), budget.read_bytes())
        # Initial coarse-only state can also acquire scheduled fine patches.
        case = source.replace('patch_count=2,', 'patch_count=0,', 1)
        run('initial_coarse_full.nml', case)
        saved = output.read_bytes(), budget.read_bytes()
        run('initial_coarse_stop.nml', options(case, "checkpoint_file='initial_coarse.chk', checkpoint_interval_steps=2, stop_after_checkpoint=.true.,"))
        run('initial_coarse_resume.nml', options(case, "restart_file='initial_coarse.chk',"))
        assert saved == (output.read_bytes(), budget.read_bytes())
    print('PASS: scheduled EB remap, budget, event restart and policy rejection')


if __name__ == '__main__':
    main()
