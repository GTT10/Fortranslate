"""Explicit CVODE SDC input, fixed/selected agreement and capability rejection."""
import argparse
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('fixed', 'selected', 'source', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--unavailable', type=Path)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='sdc_cvode_', dir=args.work.resolve()))
    text = (args.source / 'cases/reactive_boundaries_3d/sdc_reacting.nml').read_text()
    env = dict(os.environ)
    for key in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY'):
        env.pop(key, None)

    def run(label, exe, content, error=None):
        directory = work / label
        directory.mkdir()
        inp = directory / 'input.nml'
        inp.write_text(content)
        result = subprocess.run([str(exe.resolve()), str(inp)], cwd=directory, env=env,
                                capture_output=True, text=True, timeout=180)
        log = result.stdout + result.stderr
        (directory / 'run.log').write_text(log)
        output = directory / 'sdc_reacting.csv'
        if error:
            assert result.returncode != 0 and error in log, log
            assert not output.exists()
            return None
        assert result.returncode == 0, log
        residual = float(re.search(r'Maximum normalized Euler conserved change:\s+(\S+)', log).group(1))
        assert residual < 3e-12, residual
        return output, log, residual

    cvode = text.replace("time_integrator='sdc'", "time_integrator='sdc',sdc_reactor='cvode'")
    fixed, log, residual = run('fixed_cvode', args.fixed, cvode)
    assert 'SDC reactor: cvode' in log
    selected = cvode.replace("thermo_model='full_h2o2'", "thermo_model='selected'")
    selected = re.sub(r'x_h2=.*?x_n2=\.55643',
                      "composition_count=6,composition_species='H2','H','O','O2','OH','N2',"
                      'composition_mole_fractions=.29570,1e-5,1e-5,.14784,1e-5,.55643', selected)
    selected_out, selected_log, _ = run('selected_cvode', args.selected, selected)
    assert 'SDC reactor: cvode' in selected_log
    assert fixed.read_bytes() == selected_out.read_bytes(), 'Fixed/selected CVODE CSVs differ'
    native, native_log, _ = run('native_default', args.fixed, text)
    assert 'SDC reactor: native' in native_log
    with fixed.open() as stream:
        actual = list(csv.DictReader(stream))
    with native.open() as stream:
        expected = list(csv.DictReader(stream))
    for rows in (actual, expected):
        heat = [float(row['heat_release_W_m3']) for row in rows]
        assert heat and all(math.isfinite(value) for value in heat)
        assert max(abs(value) for value in heat) > 0, 'Accepted-stage reaction heat must be active'
    inert, _, _ = run('native_nonreacting', args.fixed,
                      text.replace('chemistry_enabled=.true.', 'chemistry_enabled=.false.'))
    with inert.open() as stream:
        assert all(float(row['heat_release_W_m3']) == 0 for row in csv.DictReader(stream))
    temperature_delta = max(abs(float(a['temperature']) - float(b['temperature']))
                            for a, b in zip(actual, expected, strict=True))
    assert temperature_delta < 1e-4, temperature_delta
    assert max(float(row['temperature']) for row in actual) > 1100.001, 'Reaction must heat the gas'
    run('unknown', args.fixed, cvode.replace("'cvode'", "'unknown'"), 'Unknown SDC reactor')
    run('nonreacting', args.fixed, cvode.replace('chemistry_enabled=.true.', 'chemistry_enabled=.false.'),
        'SDC CVODE requires reacting SDC')
    if args.unavailable:
        run('unavailable', args.unavailable, cvode, 'SDC CVODE backend is unavailable')
    print(f'PASS: exact fixed/selected CVODE CSVs; native temperature difference {temperature_delta:.17g} K; '
          f'conserved residual {residual:.17g}; invalid selections rejected')
    print(work)


if __name__ == '__main__':
    main()
