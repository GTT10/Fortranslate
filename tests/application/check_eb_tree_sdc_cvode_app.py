"""Check serial/MPI EB tree CVODE SDC and rank-neutral continuation."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('serial', 'mpi', 'source', 'work', 'mpiexec'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--recursive', action='store_true')
    parser.add_argument('--selected', action='store_true')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='eb_tree_cvode_', dir=args.work.resolve()))
    template = (args.source / 'cases/reactive_eb_tree_3d/sdc.nml').read_text()
    template = template.replace('time_integrator="sdc"', 'time_integrator="sdc",sdc_reactor="cvode"')
    template = template.replace('child_count=4', 'child_count=1').replace('initial_pressure=300000', 'initial_pressure=135000')
    if args.recursive:
        template = template.replace('&eb_tree', '&eb_tree subcycling=.true.,', 1)
    if args.selected:
        template = template.replace('&eb_tree', "&eb_tree composition_count=3,composition_species='H2','O2','N2',"
                                    'composition_mole_fractions=2.,1.,3.76,', 1)
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1', HWLOC_COMPONENTS='-x86',
               OMPI_MCA_hwloc_base_binding_policy='none')
    for key in ('DISPLAY', 'LD_LIBRARY_PATH', 'LIBRARY_PATH'):
        env.pop(key, None)

    def run(label, ranks=0, controls='', reactor='cvode', error=None):
        inp = work / (label + '.nml')
        output = work / (label + '.csv')
        inp.write_text(template.replace('eb_tree_sdc.csv', str(output))
                       .replace('sdc_reactor="cvode"', f'sdc_reactor="{reactor}"')
                       .replace('\n/', '\n' + controls + '\n/'))
        command = [str((args.mpi if ranks else args.serial).resolve()), str(inp)]
        if ranks:
            command = [str(args.mpiexec.resolve()), '-n', str(ranks), *command]
        result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True, timeout=180)
        log = result.stdout + result.stderr
        (work / (label + '.log')).write_text(log)
        if error:
            assert result.returncode != 0 and error in log, log
            assert not output.exists()
            return None
        assert result.returncode == 0, log
        if output.exists():
            assert 'SDC reactor: cvode' in log
            residual = float(re.search(r'Composite inventory error:\s+(\S+)', log).group(1))
            assert residual < 3e-12, residual
            return output.read_bytes()
        return None

    baseline = run('serial')
    assert run('mpi', 2) == baseline, 'Serial/MPI gas CSV differs'
    run('stop', 2, "checkpoint_file='tree.chk',checkpoint_interval_steps=1,stop_after_checkpoint=.true.,")
    assert (work / 'tree.chk').exists()
    restart = "restart_file='tree.chk',"
    assert run('restart', 4, restart) == baseline, '2-to-4-rank continuation differs'
    assert run('serial_restart', 0, restart) == baseline, 'MPI-to-serial continuation differs'
    run('changed_mpi', 4, restart, 'native', 'Invalid distributed tree checkpoint record')
    run('changed_serial', 0, restart, 'native', 'Invalid tree checkpoint record')
    print('PASS: EB tree CVODE; exact serial/MPI and rank-changing continuation; changed reactor rejected')
    print(work)


if __name__ == '__main__':
    main()
