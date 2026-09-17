#!/usr/bin/env python3
"""Check fixed/selected extruded EB full-split execution and geometry-bound restart."""
import argparse
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('fixed', 'selected', 'case-dir', 'workdir'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    template = (args.case_dir / 'reference.nml').read_text()
    section_bytes = (args.case_dir / 'section.dat').read_bytes()
    reference = None
    for kind in ('fixed', 'selected'):
        folder = (args.workdir / kind).resolve()
        folder.mkdir(parents=True, exist_ok=True)
        section = folder / 'section.dat'
        section.write_bytes(section_bytes)
        config = template
        if kind == 'selected':
            config = config.replace("thermo_model='elementary'", "thermo_model='selected'\n"
                                    " composition_count=7\n"
                                    " composition_species='H2','H','O','O2','OH','H2O','N2'\n"
                                    " composition_mole_fractions=.29570,.00001,.00001,.14784,.00001,0.,.55643")
        output = folder / 'extruded.csv'

        def run(label, content, error=None):
            inputs = folder / 'input.nml'
            inputs.write_text(content)
            result = subprocess.run([str(getattr(args, kind).resolve()), str(inputs)], cwd=folder,
                                    capture_output=True, text=True, timeout=120)
            log = result.stdout + result.stderr
            (folder / (label + '.log')).write_text(log)
            if error is None:
                assert result.returncode == 0, log
            else:
                assert result.returncode != 0 and error in log, log
            return log

        log = run('continuous', config)
        steps = int(re.search(r'Completed steps:\s+(\d+)', log)[1])
        assert steps > 1, log
        assert 'Geometry: extruded' in log and 'Molecular transport:  T' in log, log
        baseline = output.read_bytes()
        if reference is None:
            reference = baseline
        else:
            assert baseline == reference, 'Fixed/selected extruded results differ'
        save = config.replace(" output_file=", " checkpoint_file='saved.chk', checkpoint_interval_steps=1, "
                              "stop_after_checkpoint=.true.\n output_file=")
        log = run('stop', save)
        assert 'Stopped after checkpoint' in log and 'Completed steps: 1' in log, log
        restart = config.replace(" output_file=", " restart_file='saved.chk'\n output_file=")
        # Source path is not identity: the same nodal geometry may be relocated.
        (folder / 'relocated.dat').write_bytes(section_bytes)
        log = run('restart', restart.replace("section.dat", "relocated.dat"))
        assert 'Restart continuation: complete' in log, log
        assert output.read_bytes() == baseline, 'Extruded restart mismatch'
        # A valid but different section must fail checkpoint geometry comparison.
        section.write_text('10 8\n' + '\n'.join(' '.join(f'{(i/10-.5)**2+(j/8-.5)**2-.27**2:.17e}'
                           for i in range(11)) for j in range(9)) + '\n')
        run('changed_section', restart, 'Failed to restart public 3D EB state')
        assert output.read_bytes() == baseline, 'Rejected geometry changed output'
        section.write_text('9 8\n')
        run('bad_shape', config, 'Failed to build public 3D EB geometry')
        assert output.read_bytes() == baseline
        section.write_bytes(section_bytes)
        run('alias', config.replace("output_file='extruded.csv'", "output_file='section.dat'"),
            'section file must differ')
        assert section.read_bytes() == section_bytes
    print('PASS: extruded fixed/selected full split, exact restart, relocated geometry, mismatch and input guards')


if __name__ == '__main__':
    main()
