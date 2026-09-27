"""HIT Reynolds/Prandtl coefficient construction and explicit-input equivalence."""
import argparse
import json
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
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='hit_scales_', dir=args.work.resolve()))
    velocity = work / 'velocity.csv'
    rows = [(x, y, z, y, z, x) for z in (.25, .75) for y in (.25, .75) for x in (.25, .75)]
    velocity.write_text('# x,y,z,u,v,w\n' + '\n'.join(','.join(map(str, r)) for r in rows) + '\n')
    text = (args.source / 'cases/reactive_hit_3d/hit.nml').read_text()
    text = text.replace('nx=16,ny=16,nz=16', 'nx=4,ny=4,nz=4').replace('6.283185307179586', '1')
    text = text.replace('final_time=1e-6', 'final_time=1e-8')
    text = text.replace("input_file='hit_ic_4_32.dat',inres=32", f"input_file='{velocity}',inres=2")
    text = text.replace('transport_enabled=.true.',
                        'transport_enabled=.true.,species_diffusion_enabled=.false.,barodiffusion_enabled=.false.')
    controls = "&reactive_constant_transport enabled=.true.,mode='hit',lambda0=.25,reynolds_lambda0=80,prandtl=.72 /\n"
    env = dict(os.environ)
    for key in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY'):
        env.pop(key, None)

    def run(label, content, selected=False, error=None):
        directory = work / label
        directory.mkdir()
        if selected:
            content = content.replace("thermo_model='full_h2o2'", "thermo_model='selected'")
            content = re.sub(r'x_h2=0,.*?x_n2=1',
                             "composition_count=1,composition_species='N2',composition_mole_fractions=1", content)
        inp = directory / 'input.nml'
        inp.write_text(content)
        exe = args.selected if selected else args.fixed
        result = subprocess.run([str(exe.resolve()), str(inp)], cwd=directory, env=env,
                                capture_output=True, text=True, timeout=120)
        log = result.stdout + result.stderr
        (directory / 'run.log').write_text(log)
        output = directory / 'hit.csv'
        if error:
            assert result.returncode != 0 and error in log, log
            assert not output.exists()
            return
        assert result.returncode == 0, log
        residual = float(re.search(r'Maximum normalized Euler conserved change:\s+(\S+)', log).group(1))
        assert residual < 3e-12, residual
        coefficients = tuple(map(float, re.search(r'Constant viscosity/conductivity:\s+(\S+)\s+(\S+)', log).groups()))
        return output, coefficients, residual

    fixed, actual, residual = run('derived', text + controls)
    nitrogen = next(t for t in json.loads((args.source / 'mechanisms/h2o2_full.json').read_text())['thermo']
                    if t['name'] == 'N2')
    gas_r = 8314.46261815324 / nitrogen['molecular_weight']
    cp = gas_r * sum(a * 400**i for i, a in enumerate(nitrogen['low_coefficients'][:5]))
    rho = 101325 / (gas_r * 400)
    rms = .1 * math.sqrt(cp / (cp - gas_r) * gas_r * 400) / math.sqrt(3)
    mu = rho * rms * .25 / 80
    expected = mu, mu * cp / .72
    error = max(abs(a - b) / b for a, b in zip(actual, expected))
    assert error < 3e-14, (actual, expected)
    explicit = f'&reactive_constant_transport enabled=.true.,viscosity={actual[0]:.17g},conductivity={actual[1]:.17g} /\n'
    equivalent, _, _ = run('specified', text + explicit)
    selected, _, _ = run('selected', text + controls, selected=True)
    assert fixed.read_bytes() == equivalent.read_bytes() == selected.read_bytes()
    _, half, _ = run('double_reynolds', text + controls.replace('reynolds_lambda0=80', 'reynolds_lambda0=160'))
    assert half == tuple(v / 2 for v in actual)
    _, normalized, _ = run('different_file_normalization', text.replace('uin_norm=1', 'uin_norm=2') + controls)
    assert normalized == actual, 'File normalization must not change target Reynolds coefficients'
    run('zero_reynolds', text + controls.replace('reynolds_lambda0=80', 'reynolds_lambda0=0'),
        error='HIT transport parameters must be positive')
    run('mixed_controls', text + controls.replace("mode='hit'", "mode='hit',viscosity=.001"),
        error='HIT mode cannot also specify transport coefficients')
    run('missing_velocity', re.sub(r'&reactive_hit\b.*?^/\s*$', '', text, flags=re.S | re.M) + controls,
        error='HIT transport requires positive initialized HIT scales')
    run('overflow', text + controls.replace('lambda0=.25', 'lambda0=1e308'),
        error='HIT transport coefficient overflow')
    print(f'PASS: coefficient relative error {error:.17g}; specified/selected exact CSV; '
          f'Reynolds scaling and rejection; residual {residual:.17g}')
    print(work)


if __name__ == '__main__':
    main()
