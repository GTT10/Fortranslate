#!/usr/bin/env python3
"""Selected/fixed elementary tree parity through prescribed regridding."""
import argparse
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--selected', type=Path, required=True)
    parser.add_argument('--fixed', type=Path, required=True)
    parser.add_argument('--case-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    work = args.work_dir.resolve()
    case = args.case_dir.resolve()
    for name in ['fixed', 'selected']:
        (work / name).mkdir(parents=True, exist_ok=True)
    config = (case / 'hotspot.nml').read_text().replace("thermo_model='elementary'", "thermo_model='selected'")
    config = config.replace("  output_file=", "  composition_count=7\n"
                            "  composition_species='H2','H','O','O2','OH','H2O','N2'\n"
                            "  composition_mole_fractions=.29570,.00001,.00001,.14784,.00001,0.0,.55643\n"
                            "  output_file=")
    selected_input = work / 'selected.nml'
    selected_input.write_text(config)

    def run(binary, inputs, directory, label, success=True):
        result = subprocess.run([str(binary.resolve()), *map(str, inputs)], cwd=directory,
                                text=True, capture_output=True, timeout=180)
        log = result.stdout + result.stderr
        (directory / f'{label}.log').write_text(log)
        assert (result.returncode == 0) == success, log
        return log

    tail = [case / 'branched.dat', case / 'regrid.dat']
    run(args.fixed, [case / 'hotspot.nml', *tail], work / 'fixed', 'run')
    log = run(args.selected, [selected_input, *tail], work / 'selected', 'run')
    assert 'Bundle SHA-256:' in log and 'Chemistry integrator:' in log
    for name in ['tree_hotspot.csv', 'tree_hotspot.csv.budget.csv']:
        assert (work / 'fixed' / name).read_bytes() == (work / 'selected' / name).read_bytes(), name
    baseline = (work / 'selected' / 'tree_hotspot.csv').read_bytes()
    invalid = work / 'invalid.nml'
    invalid.write_text(config.replace("'H2','H'", "'unknown','H'"))
    run(args.selected, [invalid, *tail], work / 'selected', 'invalid_composition', False)
    assert (work / 'selected' / 'tree_hotspot.csv').read_bytes() == baseline
    invalid.write_text(config.replace('initial_temperature=1100.0', 'initial_temperature=7000.0'))
    run(args.selected, [invalid, *tail], work / 'selected', 'invalid_temperature', False)
    assert (work / 'selected' / 'tree_hotspot.csv').read_bytes() == baseline
    print('PASS: selected/fixed exact tree CSV and budgets through regridding; startup rejections')


if __name__ == '__main__':
    main()
