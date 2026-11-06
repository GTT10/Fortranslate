#!/usr/bin/env python3
"""Fresh reacting MOL initialization and rank-neutral application continuation."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('serial', 'mpi', 'mpiexec', 'work'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--transport', action='store_true', help='Enable application molecular transport')
    parser.add_argument('--subcycling', type=int, choices=(0, 1), default=1,
                        help='Common (0) or recursive (1, default) level timesteps')
    args = parser.parse_args()
    args.serial = args.serial.resolve(strict=True)
    args.mpi = args.mpi.resolve(strict=True)
    args.mpiexec = args.mpiexec.resolve(strict=True)
    root = Path(__file__).resolve().parents[1]
    case = root / 'cases/amr_tree_reactive_3d'
    template = (case / 'mol_initial_reactions.nml').read_text()
    if args.transport:
        assert template.count('transport_enabled=.false.') == 1
        template = template.replace('transport_enabled=.false.', 'transport_enabled=.true.')
    if args.subcycling == 0:
        assert template.count('subcycling=.true.') == 1
        template = template.replace('subcycling=.true.', 'subcycling=.false.')
    args.work.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='initial_mol_', dir=args.work.resolve()))
    shutil.copyfile(case / 'mol_initial_reactions.dat', work / 'topology.dat')
    print(work, flush=True)
    report = {'status': 'running', 'commands': [], 'checker_sha256': digest(Path(__file__)),
              'case_sha256': digest(case / 'mol_initial_reactions.nml'),
              'topology_sha256': digest(work / 'topology.dat'),
              'serial_sha256': digest(args.serial), 'mpi_sha256': digest(args.mpi)}
    report['transport_enabled'] = args.transport
    report['subcycling'] = args.subcycling
    env = dict(os.environ, HWLOC_COMPONENTS='-gl,-x86', OMP_NUM_THREADS='1',
               OMPI_MCA_rmaps_base_oversubscribe='1')
    marker = 'Initialized MOL reaction history for first accepted step'

    def save():
        (work / 'checks.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')

    def run(label, ranks, extra='', restart=False):
        config = work / (label + '.nml')
        config.write_text(template + extra)
        command = [str(args.serial if ranks == 0 else args.mpi), str(config), str(work / 'topology.dat')]
        if ranks:
            command = [str(args.mpiexec), '-n', str(ranks), *command]
        entry = {'label': label, 'command': command, 'cwd': str(work), 'input_sha256': digest(config)}
        report['commands'].append(entry)
        save()
        start = time.monotonic()
        try:
            with (work / (label + '.log')).open('w') as stream:
                result = subprocess.run(command, cwd=work, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=120)
            entry['returncode'] = result.returncode
            assert result.returncode == 0, label
            log = (work / (label + '.log')).read_text()
            assert log.count(marker) == (0 if restart else 1), (label, 'initialization lifecycle')
            if label != 'stop_mpi2':
                entry['output_sha256'] = {}
                for name in ('result.csv', 'result.csv.budget.csv'):
                    output = work / name
                    snapshot = work / (label + '_' + name)
                    shutil.copyfile(output, snapshot)
                    entry['output_sha256'][name] = digest(snapshot)
            return log
        finally:
            entry['seconds'] = time.monotonic() - start
            save()

    outputs = [work / 'result.csv', work / 'result.csv.budget.csv']
    try:
        log = run('serial', 0)
        assert 'Steps: 2' in log
        baseline = [p.read_bytes() for p in outputs]
        for p in outputs:
            shutil.copyfile(p, work / ('serial_' + p.name))
        with outputs[1].open() as stream:
            budget = list(csv.DictReader(stream))
        with outputs[0].open() as stream:
            rows = list(csv.DictReader(stream))
        assert rows and all(math.isfinite(float(v)) for r in rows for v in r.values())
        assert len(rows) == 71
        assert all(abs(float(r['time']) - 2.e-7) < 1.e-19 for r in rows)
        mechanism = json.loads((root / 'mechanisms/lidryer.json').read_text())
        names = mechanism['species']
        fields = ['rho', 'rhou', 'rhov', 'rhow', 'rhoE'] + ['rhoY_' + n for n in names]
        volumes = [math.prod(float(r[f]) for f in ('dx', 'dy', 'dz')) for r in rows]
        assert len(budget) == 14
        residuals = []
        for field, row in zip(fields, budget):
            initial, final = float(row['initial']), float(row['final'])
            changes = sum(float(row[k]) for k in ('boundary_change', 'source_change', 'energy_correction'))
            inventory = math.fsum(float(r[field]) * v for r, v in zip(rows, volumes))
            # Signed momentum can sum to zero despite finite local motion.
            # Normalize its cancellation residual by the absolute inventory.
            absolute_inventory = math.fsum(abs(float(r[field])) * v for r, v in zip(rows, volumes))
            scale = max(abs(initial), abs(final), abs(changes), absolute_inventory, 1.e-30)
            assert abs(inventory - final) <= 1.e-9 * scale, (field, 'CSV inventory')
            # source_change is the external-source ledger; chemical species
            # conversion is checked through elements, not individual species.
            if not field.startswith('rhoY_'):
                assert abs(final - initial - changes) <= 1.e-9 * scale, row
                residuals.append(abs(final - initial - changes) / scale)
        report['maximum_scaled_budget_residual'] = max(residuals)
        report['element_relative_drift'] = {}
        for element in ('H', 'O', 'N'):
            inventories = {}
            for label in ('initial', 'final'):
                inventories[label] = math.fsum(float(row[label]) * sp['composition'].get(element, 0.) /
                                              sp['molecular_weight'] for row, sp in zip(budget[5:], mechanism['thermo']))
            drift = abs(inventories['final'] / inventories['initial'] - 1)
            assert drift <= 1.e-9, (element, drift)
            report['element_relative_drift'][element] = drift
        assert all(float(r['rhoY_' + n]) >= 0 for r in rows for n in names)
        assert all(abs(sum(float(r['rhoY_' + n]) for n in names) / float(r['rho']) - 1) <= 1.e-12 for r in rows)
        consumed = 1 - float(budget[5]['final']) / float(budget[5]['initial'])
        assert consumed > 1.e-6, 'resolved chemical change required'
        report['hydrogen_consumed_fraction'] = consumed
        run('mpi2', 2)
        assert [p.read_bytes() for p in outputs] == baseline, 'serial/MPI final fields and budgets'
        stop = "\n&amr_tree_checkpoint checkpoint_file='state.chk',checkpoint_interval_steps=1,stop_after_checkpoint=.true. /\n"
        log = run('stop_mpi2', 2, stop)
        assert 'Stopped after checkpoint at step 1' in log
        assert [p.read_bytes() for p in outputs] == baseline, 'stop must not replace final output'
        checkpoint = digest(work / 'state.chk')
        resume = "\n&amr_tree_checkpoint restart_file='state.chk' /\n"
        for label, ranks in (('resume_mpi3', 3), ('resume_serial', 0)):
            log = run(label, ranks, resume, restart=True)
            assert 'Steps: 2' in log
            assert digest(work / 'state.chk') == checkpoint
            assert [p.read_bytes() for p in outputs] == baseline, label
        report['checkpoint_sha256'] = checkpoint
        report['status'] = 'passed'
        print('PASS: fresh MOL history, serial/MPI exact output, 2-to-3-rank and serial restart without reseeding', flush=True)
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = repr(error)
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
