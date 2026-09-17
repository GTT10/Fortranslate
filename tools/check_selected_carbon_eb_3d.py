#!/usr/bin/env python3
"""Check the small single-level carbon EB continuation fixture."""
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
    exe, work = args.exe.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    source = args.case.read_text()

    def run(name, content, success=True):
        (work / name).write_text(content)
        result = subprocess.run([str(exe), name], cwd=work, text=True, capture_output=True, timeout=90)
        log = result.stdout + result.stderr
        (work / (name + '.log')).write_text(log)
        assert (result.returncode == 0) == success, log
        return log

    log = run('full.nml', source)
    assert 'Operator sequence: R-T-H-T-R' in log
    steps = int(next(line.split(':')[1] for line in log.splitlines() if line.startswith('Completed steps:')))
    assert steps > 1, 'fixture must exercise continuation before its final time'
    output = work / 'methanol_eb.csv'
    saved = output.read_bytes()
    rows = list(csv.DictReader(saved.decode().splitlines()))
    assert len(rows) == 8
    assert len([key for key in rows[0] if key.startswith('rhoY_')]) == 21
    for row in rows:
        assert all(math.isfinite(float(value)) for value in row.values())
        assert float(row['rho']) > 0 and float(row['temperature']) > 0
        assert math.isclose(sum(float(v) for k, v in row.items() if k.startswith('rhoY_')),
                            float(row['rho']), rel_tol=5e-11)
    stopped = source.replace('\n/', "\n checkpoint_file='carbon.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true.,\n/")
    run('stopped.nml', stopped)
    resume = source.replace('\n/', "\n restart_file='carbon.chk',\n/")
    run('resume.nml', resume)
    assert output.read_bytes() == saved
    resumed_log = (work / 'resume.nml.log').read_text()
    budget_lines = lambda text: [line for line in text.splitlines() if line.startswith('EB inventory')]
    assert budget_lines(log) == budget_lines(resumed_log), 'cumulative restart budgets differ'
    checkpoint = work / 'carbon.chk'
    original = checkpoint.read_text()
    assert original.splitlines()[1].split()[0] == '6'
    lines = original.splitlines()
    index = lines.index('ELEMENT_CONTEXT')
    assert int(lines[index + 1]) == 5
    assert lines[index + 2] == 'C'
    lines[index + 2] = 'WRONG_C'
    checkpoint.write_text('\n'.join(lines) + '\n')
    run('bad_element.nml', resume, False)
    assert output.read_bytes() == saved
    checkpoint.write_text(original)
    lines = original.splitlines()
    exchange = lines.index('EXCHANGES')
    values = lines[exchange + 1].split()
    values[0] = '1.0e-4'
    lines[exchange + 1] = ' '.join(values)
    corrupted = '\n'.join(lines) + '\n'
    checkpoint.write_text(corrupted)
    run('bad_budget.nml', resume, False)
    assert output.read_bytes() == saved
    assert checkpoint.read_text() == corrupted
    checkpoint.write_text(original)
    print('PASS: single-level carbon EB state, elemental checkpoint and exact continuation')


if __name__ == '__main__':
    main()
