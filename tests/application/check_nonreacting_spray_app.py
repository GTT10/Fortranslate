"""Check reaction-free NC7H16 spray integration, mass transfer and input guards."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'source', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='nonreacting_', dir=args.work.resolve()))
    base = (args.source / 'cases/spray_3d/heptane_nonreacting.nml').read_text()
    env = dict(os.environ)
    for key in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY'):
        env.pop(key, None)
    for name, text, error in (
        ('evaporation', base, None),
        ('chemistry', base.replace('chemistry_enabled=.false.', 'chemistry_enabled=.true.'),
         'Nonreacting spray mechanism requires chemistry disabled'),
        ('cold', base.replace('gas_temperature=400', 'gas_temperature=273'),
         'Initial gas state is outside the mechanism thermodynamic range'),
    ):
        directory = work / name
        directory.mkdir()
        (directory / 'input.nml').write_text(text)
        result = subprocess.run([str(args.exe.resolve()), 'input.nml'], cwd=directory, env=env,
                                capture_output=True, text=True, timeout=120)
        log = result.stdout + result.stderr
        (directory / 'run.log').write_text(log)
        if error:
            assert result.returncode != 0 and error in log, log
            assert not (directory / 'heptane.gas.csv').exists()
        else:
            assert result.returncode == 0, log
    directory = work / 'evaporation'
    with (directory / 'heptane.history.csv').open() as stream:
        history = list(csv.DictReader(stream))
    assert len(history) == 3 and int(history[-1]['steps']) == 2
    initial, final = history[0], history[-1]
    loss = float(initial['liquid_mass']) - float(final['liquid_mass'])
    assert loss > 1e-16, loss
    with (directory / 'heptane.gas.csv').open() as stream:
        vapor = sum(float(row['rhoY_NC7H16']) for row in csv.DictReader(stream)) * (.004/4)**3
    assert abs(vapor-loss)/loss < 1e-9, (vapor, loss)
    mass_error = max(abs(float(row['mass_error'])) for row in history) / float(initial['system_mass'])
    energy_error = max(abs(float(row['energy_error'])) for row in history) / abs(float(initial['total_energy']))
    assert mass_error < 1e-12 and energy_error < 1e-12, (mass_error, energy_error)
    report = dict(evaporated_mass=loss, vapor_mass=vapor, normalized_mass_error=mass_error,
                  normalized_energy_error=energy_error)
    (work / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS:', report)
    print(work)


if __name__ == '__main__':
    main()
