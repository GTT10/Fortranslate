#!/usr/bin/env python3
"""Exercise selected two-level EB composition, inventories and restart identity."""
import argparse
import csv
import hashlib
import math
import re
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    for name in ('exe', 'case', 'bundle', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--species-count', type=int, default=10)
    parser.add_argument('--fuel')
    args = parser.parse_args()
    exe, work = args.exe.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    source = args.case.read_text()
    output = work / 'eb_selected.csv'
    budget = work / 'eb_selected.csv.budget.csv'

    def run(name, content, success=True):
        (work / name).write_text(content)
        result = subprocess.run([str(exe), name], cwd=work, capture_output=True, text=True, timeout=90)
        (work / (name + '.log')).write_text(result.stdout + result.stderr)
        assert (result.returncode == 0) == success, result.stdout + result.stderr

    def options(content, extra):
        return content.replace('\n/', '\n ' + extra + '\n/')

    for count, subcycle in ((2, False), (1, True), (2, True)):
        case = source.replace('patch_count=2', f'patch_count={count}, subcycling={".true." if subcycle else ".false."}')
        count = f'{count}_{int(subcycle)}'
        run(f'full{count}.nml', case)
        saved = output.read_bytes(), budget.read_bytes()
        rows = list(csv.DictReader(output.read_text().splitlines()))
        assert len(rows) == 104
        fields = [key for key in rows[0] if key == 'rho' or key.startswith(('rhou', 'rhov', 'rhow', 'rhoE', 'rhoY_'))]
        assert len([key for key in fields if key.startswith('rhoY_')]) == args.species_count
        totals = {key: 0.0 for key in fields}
        for row in rows:
            assert all(math.isfinite(float(value)) for value in row.values())
            assert float(row['rho']) > 0 and float(row['temperature']) > 0
            species = sum(float(row[key]) for key in fields if key.startswith('rhoY_'))
            assert math.isclose(species, float(row['rho']), rel_tol=1e-10)
            volume = float(row['kappa']) / (4 if int(row['level']) == 0 else 8) ** 3
            for key in fields:
                totals[key] += volume * float(row[key])
        budgets = list(csv.DictReader(budget.read_text().splitlines()))
        assert len(budgets) == len(fields)
        for key, entry in zip(fields, budgets):
            assert math.isclose(totals[key], float(entry['final']), rel_tol=3e-12, abs_tol=3e-12)
            residual = float(entry['final']) - float(entry['initial']) - sum(float(entry[k]) for k in ('boundary', 'wall', 'chemistry'))
            assert abs(residual) <= 5e-10 * max(1.0, abs(float(entry['initial'])), abs(float(entry['final'])))
        assert any(abs(float(entry['chemistry'])) > 1e-15 for entry in budgets[5:])
        if args.fuel:
            fuel_index = fields.index('rhoY_' + args.fuel)
            assert abs(float(budgets[fuel_index]['chemistry'])) > 1e-15
        run(f'stop{count}.nml', options(case, "checkpoint_file='selected.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true.,"))
        checkpoint = work / 'selected.chk'
        original = checkpoint.read_text()
        sha = hashlib.sha256(args.bundle.read_bytes()).hexdigest()
        assert sha in original
        resumed = options(case, "restart_file='selected.chk',")
        run(f'resume{count}.nml', resumed)
        assert saved == (output.read_bytes(), budget.read_bytes())
        wrong_mode = resumed.replace('subcycling=.true.' if subcycle else 'subcycling=.false.',
                                     'subcycling=.false.' if subcycle else 'subcycling=.true.')
        run(f'wrong_mode{count}.nml', wrong_mode, False)
        assert saved == (output.read_bytes(), budget.read_bytes())
        checkpoint.write_text(original.replace(sha, '0' * 64))
        run(f'wrong_sha{count}.nml', resumed, False)
        assert saved == (output.read_bytes(), budget.read_bytes())
        metadata = original.index('ELEMENT_CONTEXT\n')
        lines = original[metadata:].splitlines()
        assert int(lines[1]) >= 3
        lines[2] = 'WRONG_ELEMENT'
        checkpoint.write_text(original[:metadata] + '\n'.join(lines) + '\n')
        run(f'wrong_element{count}.nml', resumed, False)
        assert saved == (output.read_bytes(), budget.read_bytes())
        lines = original[metadata:].splitlines()
        atoms = lines[3].split()
        atoms[0] = str(float(atoms[0]) + 1)
        lines[3] = ' '.join(atoms)
        checkpoint.write_text(original[:metadata] + '\n'.join(lines) + '\n')
        run(f'wrong_atoms{count}.nml', resumed, False)
        assert saved == (output.read_bytes(), budget.read_bytes())
        checkpoint.write_text(original)
        changed_mix = re.sub(r'composition_mole_fractions=[^\n]+',
                             'composition_mole_fractions=0.30,0.15,0.55,', resumed)
        run(f'wrong_mix{count}.nml', changed_mix, False)
        assert saved == (output.read_bytes(), budget.read_bytes())
    run('bad_species.nml', re.sub(r"composition_species='[^']+'", "composition_species='BAD'", source), False)
    print('PASS: selected EB single/multiple patches, chemistry, inventories and bound restart')


if __name__ == '__main__':
    main()
