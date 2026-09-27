"""Check cut-cell CVODE SDC and reactor-specific checkpoint continuation."""
import argparse
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
    work = Path(tempfile.mkdtemp(prefix='eb_cvode_', dir=args.work.resolve()))
    template = (args.source / 'cases/reactive_eb_3d/sdc.nml').read_text()
    env = dict(os.environ)
    for key in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY'):
        env.pop(key, None)

    for selected in (False, True):
        folder = work / ('selected' if selected else 'fixed')
        folder.mkdir()
        text = template
        if selected:
            text = text.replace('&reactive_eb_3d', "&reactive_eb_3d thermo_model='selected', "
                                "composition_count=3,composition_species='H2','O2','N2',"
                                'composition_mole_fractions=2.,1.,3.76,')
        cvode = text.replace('time_integrator="sdc"', 'time_integrator="sdc",sdc_reactor="cvode"')

        def run(label, content, controls='', error=None):
            inp = folder / (label + '.nml')
            output = folder / (label + '.csv')
            content = content.replace('reactive_eb_3d_sdc.csv', str(output))
            inp.write_text(content.replace('\n/', '\n' + controls + '\n/'))
            result = subprocess.run([str((args.selected if selected else args.fixed).resolve()), str(inp)],
                                    cwd=folder, env=env, capture_output=True, text=True, timeout=180)
            log = result.stdout + result.stderr
            (folder / (label + '.log')).write_text(log)
            if error:
                assert result.returncode != 0 and error in log, log
                assert not output.exists(), 'Rejected input/restart published output'
                return None, log
            assert result.returncode == 0, log
            return output.read_bytes(), log

        baseline, log = run('continuous', cvode)
        assert 'SDC reactor: cvode' in log
        inventory = re.findall(r'^EB inventory initial/final/boundary/wall/chemical/forcing/residual:\s+(.+)$',
                               log, re.MULTILINE)
        assert inventory, log
        values = [[float(x) for x in row.split()[1:]] for row in inventory]
        mass_scale = abs(values[0][0])
        residual = max(abs(row[-1]) / max(mass_scale, abs(row[0])) for row in values)
        assert residual < 3e-12, residual
        assert any(abs(row[4]) > 0 for row in values[5:]), 'Chemistry must convert species'
        _, stopped = run('stop', cvode, "checkpoint_file='cvode.chk',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,")
        assert 'Stopped after checkpoint' in stopped
        assert 'SDC_CVODE_SOURCES' in (folder / 'cvode.chk').read_text()
        resumed, resumed_log = run('resumed', cvode, "restart_file='cvode.chk',")
        assert resumed == baseline, 'CVODE checkpoint continuation differs'
        assert 'Restarted from checkpoint:' in resumed_log
        run('changed_native', text, "restart_file='cvode.chk',", 'Failed to restart public 3D EB state')
        run('native_stop', text, "checkpoint_file='native.chk',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,")
        assert 'SDC_SOURCES' in (folder / 'native.chk').read_text()
        run('changed_cvode', cvode, "restart_file='native.chk',", 'Failed to restart public 3D EB state')
        run('unknown', cvode.replace('sdc_reactor="cvode"', 'sdc_reactor="unknown"'), error='Unknown EB SDC reactor')
        print(f'PASS: {folder.name}: active CVODE, exact restart, both reactor changes rejected; residual {residual:.17g}')
    print(work)


if __name__ == '__main__':
    main()
