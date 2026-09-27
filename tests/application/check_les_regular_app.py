"""Regular gas LES SDC input, reacting fixed/selected agreement and rejection."""
import argparse
import csv
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
    work = Path(tempfile.mkdtemp(prefix='les_regular_', dir=args.work.resolve()))
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

    text = text.replace('final_time=1e-8', 'final_time=1e-10')
    text = text.replace("boundary_condition='periodic'", "boundary_condition='periodic',boundary_faces='no_slip_wall','no_slip_wall','no_slip_wall','no_slip_wall','periodic','periodic'")
    cvode = text.replace("time_integrator='sdc'", "time_integrator='sdc',sdc_reactor='cvode'")
    disabled, _, _ = run('disabled', args.fixed, cvode)
    maximum = 0.0
    for model in ('smagorinsky', 'wale', 'vreman', 'dynamic'):
        content = cvode + f"\n&reactive_les enabled=.true.,model='{model}' /\n"
        fixed, log, residual = run(model, args.fixed, content)
        assert f'LES model: {model}' in log
        assert fixed.read_bytes() != disabled.read_bytes(), 'LES must change the gas'
        maximum = max(maximum, residual)
        selected = content.replace("thermo_model='full_h2o2'", "thermo_model='selected'")
        selected = re.sub(r'x_h2=.*?x_n2=\.55643',
                          "composition_count=6,composition_species='H2','H','O','O2','OH','N2',"
                          'composition_mole_fractions=.29570,1e-5,1e-5,.14784,1e-5,.55643', selected)
        selected_out, _, _ = run(model + '_selected', args.selected, selected)
        assert fixed.read_bytes() == selected_out.read_bytes(), 'Fixed/selected LES CSVs differ'
    run('invalid_model', args.fixed, cvode + "\n&reactive_les enabled=.true.,model='invalid' /\n",
        'Failed to compute additional diffusion timestep')
    run('split', args.fixed, text.replace("time_integrator='sdc'", "time_integrator='split'") +
        "\n&reactive_les enabled=.true. /\n", 'LES requires SDC')
    if args.unavailable:
        run('unavailable', args.unavailable, text + "\n&reactive_les enabled=.true. /\n",
            'Regular LES backend is unavailable')
    print(f'PASS: four active LES closures with CVODE; exact fixed/selected CSVs; residual {maximum:.17g}')
    print(work)


if __name__ == '__main__':
    main()
