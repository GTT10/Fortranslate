"""AMR CVODE SDC serial/MPI equality and changed-rank checkpoint continuation."""
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
    for name in ('serial', 'mpi', 'source', 'work'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--recursive', action='store_true')
    parser.add_argument('--selected', action='store_true')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='amr_cvode_', dir=args.work.resolve()))
    template = (args.source / 'cases/amr_tree_reactive_3d/sdc.nml').read_text()
    template = template.replace("time_integrator='sdc'", "time_integrator='sdc',sdc_reactor='cvode'")
    template = template.replace('cfl=.004', 'cfl=.0005')
    if args.selected:
        template = template.replace("thermo_model='full_h2o2'", "thermo_model='selected'")
        template = re.sub(r'x_h2=.*?x_n2=\.55643',
                          "composition_count=6,composition_species='H2','H','O','O2','OH','N2',"
                          'composition_mole_fractions=.29570,1e-5,1e-5,.14784,1e-5,.55643', template)
    if args.recursive:
        template += '\n&amr_tree_timestepping subcycling=.true. /\n'
    topology = work / 'tree.dat'
    topology.write_text('1\n1 2 2 2 2 2 2 2\n')
    env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe='1', HWLOC_COMPONENTS='-x86',
               OMPI_MCA_hwloc_base_binding_policy='none')
    for key in ('LD_LIBRARY_PATH', 'LIBRARY_PATH', 'DISPLAY'):
        env.pop(key, None)

    def run(label, content, ranks=0, error=None):
        inp = work / (label + '.nml')
        output = work / (label + '.csv')
        inp.write_text(content.replace('sdc_tree.csv', str(output)))
        command = [str((args.mpi if ranks else args.serial).resolve()), str(inp), str(topology)]
        if ranks:
            command = [args.mpiexec, '-n', str(ranks), *command]
        result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True, timeout=180)
        log = result.stdout + result.stderr
        (work / (label + '.log')).write_text(log)
        if error:
            assert result.returncode != 0 and error in log, log
            assert not output.exists(), 'Rejected restart published output'
            return None, log
        assert result.returncode == 0, log
        outputs = None
        if output.exists():
            outputs = [output.read_bytes(), Path(str(output) + '.budget.csv').read_bytes()]
        return outputs, log

    baseline, log = run('serial', template)
    assert 'SDC reactor: cvode' in log
    assert int(re.search(r'Steps:\s+(\d+)', log).group(1)) > 1, log
    with (work / 'serial.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert max(float(row['temperature']) for row in rows) > 1100.00001, 'Chemistry must heat the gas'
    heat = [float(row['heat_release_W_m3']) for row in rows]
    assert all(math.isfinite(value) for value in heat) and any(value != 0 for value in heat)
    with (work / 'serial.csv.budget.csv').open() as stream:
        ledger = list(csv.DictReader(stream))
    residual = max(abs(float(row['final']) - float(row['initial']) - float(row['boundary_change'])
                       - float(row['source_change'])) / max(1., abs(float(row['initial'])))
                   for row in ledger[:5])
    assert residual < 3e-12, residual
    parallel, _ = run('mpi', template, 2)
    assert parallel == baseline, 'Serial/MPI output differs'
    save = "\n&amr_tree_checkpoint checkpoint_file='tree.chk', checkpoint_interval_steps=1, stop_after_checkpoint=.true. /\n"
    restart = "\n&amr_tree_checkpoint restart_file='tree.chk' /\n"
    _, log = run('stop', template + save, 2)
    assert 'Stopped after checkpoint at step 1' in log, log
    assert (work / 'tree.chk').read_bytes().startswith(b'PELEF_AMR_TREE_3D_CHECKPOINT_V5\n')
    resumed, _ = run('restart', template + restart, 4)
    assert resumed == baseline, '2-to-4-rank restart differs'
    resumed, _ = run('serial_restart', template + restart)
    assert resumed == baseline, 'MPI-to-serial restart differs'
    # Make a genuine V3 legacy archive by omitting the new scalar payload.
    checkpoint = work / 'tree.chk'
    modern = checkpoint.read_bytes()
    lines = modern.decode().splitlines()
    context_bytes = int(lines[1])
    cursor, consumed = 2, 0
    while consumed < context_bytes:
        consumed += len(lines[cursor]) + 1
        cursor += 1
    assert consumed == context_bytes
    count, nvar, nx, ny, nz = map(int, lines[cursor].split())
    shapes = [(nx, ny, nz)]
    cursor += 5  # node header, clock, ledger, forcing, reaction presence
    for _ in range(count - 1):
        parent, il, iu, jl, ju, kl, ku, ratio = map(int, lines[cursor].split())
        assert 1 <= parent <= len(shapes)
        shapes.append(((iu - il + 1) * ratio, (ju - jl + 1) * ratio, (ku - kl + 1) * ratio))
        cursor += 1
    heat_cells = sum(math.prod(shape) for shape in shapes)
    assert lines[-1] == 'END_AMR_TREE_CHECKPOINT'
    legacy = lines[:-1 - heat_cells] + [lines[-1]]
    legacy[0] = 'PELEF_AMR_TREE_3D_CHECKPOINT_V3'
    checkpoint.write_text('\n'.join(legacy) + '\n')
    resumed, _ = run('legacy_restart', template + restart, 2)
    assert resumed == baseline, 'Legacy archive cold diagnostic differs after accepted stage'
    checkpoint.write_bytes(modern)
    run('changed_reactor', template.replace("sdc_reactor='cvode'", "sdc_reactor='native'") + restart,
        4, 'Invalid MPI tree restart checkpoint/context')
    print(f'PASS: active CVODE; exact serial/MPI and changed-rank restart; reactor mismatch rejection; residual {residual:.17g}')
    print(work)


if __name__ == '__main__':
    main()
