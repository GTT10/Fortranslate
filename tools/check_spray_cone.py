#!/usr/bin/env python3
"""Check nonuniform reactive spray/LES/injection and exact coupled restart.

This tests integration and conservation, not experimental spray accuracy.
The selected executable must use the pinned FFCM1_Red reference mechanism.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile

import cantera as ct
import numpy as np


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(encoding='utf-8') as stream:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check(exe: Path, source: Path, case: Path, work: Path) -> None:
    require(ct.__version__ == '3.2.0', 'Reference requires Cantera 3.2.0')
    original = case.read_text(encoding='utf-8')
    require(original.rstrip().endswith('/'), 'Expected a spray namelist')

    def run(name: str, extra: str = '') -> None:
        # Later namelist assignments override only run paths/stop controls.
        text = original.rstrip()[:-1] + f"\n output_prefix='{work/name}',checkpoint_file='',\n{extra}\n/\n"
        inp = work/(name+'.nml')
        inp.write_text(text, encoding='utf-8')
        result = subprocess.run([str(exe), str(inp)], cwd=work, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
        (work/(name+'.log')).write_text(result.stdout, encoding='utf-8')
        require(result.returncode == 0, result.stdout)

    run('continuous')
    run('before', f"stop_after_steps=10,checkpoint_file='{work/'split.chk'}',")
    run('after', f"restart_file='{work/'split.chk'}',")
    for suffix in ('gas.csv', 'parcels.csv'):
        require((work/f'continuous.{suffix}').read_bytes() == (work/f'after.{suffix}').read_bytes(),
                'Reactive restart differs in '+suffix)
    history = read_rows(work/'continuous.history.csv')
    tail = read_rows(work/'after.history.csv')
    require([row for row in history if row['steps'] >= tail[0]['steps']] == tail,
            'Reactive restart changes the output history')
    field = read_rows(work/'continuous.gas.csv')
    parcels = read_rows(work/'continuous.parcels.csv')
    require(len(field) == 27, 'Expected 3x3x3 cells')
    require(history[-1]['time'] == 3e-6, 'Cone did not reach the requested time')
    require(history[-1]['steps'] > tail[0]['steps'] > 0, 'Restart not exercised')
    require(all(np.isfinite(v) for row in history+field+parcels for v in row.values()), 'Nonfinite output')
    injected = history[-1]['injected_mass']
    require(abs(injected-1e-10) < 1e-22, 'Incorrect scheduled injection mass')
    require(len(parcels) == 6, 'Initial and two scheduled parcel pulses were not retained')
    mass_error = max(abs(row['mass_error'])/row['system_mass'] for row in history)
    energy_error = max(abs(row['energy_error'])/max(abs(row['total_energy']),1e-12) for row in history)
    momentum_error = max(max(abs(row[k]) for k in ('px_error','py_error','pz_error')) for row in history)
    require(mass_error < 2e-11 and energy_error < 2e-11 and momentum_error < 2e-17,
            'Coupled nonuniform conservation failed')
    gas = ct.Solution(str(source), 'gas')
    gas.TPX = 2200., 5e5, {'O2':1.,'N2':3.76}
    volume = .004**3
    matrix = np.array([[gas.n_atoms(k,e)/gas.molecular_weights[k] for k in range(gas.n_species)]
                       for e in gas.element_names])
    initial_species = gas.Y*gas.density*volume
    initial_species[gas.species_index('CH3OH')] += 1e-9+injected
    actual_species = np.mean([[row['rhoY_'+n] for n in gas.species_names] for row in field],axis=0)*volume
    remaining = sum(p['mass']*p['multiplicity'] for p in parcels)
    actual_species[gas.species_index('CH3OH')] += remaining
    expected_elements = matrix @ initial_species
    element_delta = np.abs(matrix @ actual_species-expected_elements)
    present = expected_elements > 0
    element_error = float(np.max(element_delta[present]/expected_elements[present]))
    # Initially absent elements have no relative reference. The gas solver closes
    # its last species (He here); bound its roundoff using total elemental moles.
    absent_error = float(np.max(element_delta[~present], initial=0.0))
    require(element_error < 2e-8, 'Nonuniform element inventory mismatch')
    require(absent_error < 64*np.finfo(float).eps*np.sum(expected_elements),
            'An initially absent element exceeds the absolute roundoff bound')
    temperature_range = max(row['temperature'] for row in field)-min(row['temperature'] for row in field)
    require(temperature_range > 1., 'Cone test is not spatially nonuniform')
    require(remaining < 1e-9+injected, 'Cone did not evaporate')
    products = sum(actual_species[gas.species_index(n)] for n in ('CO','CO2'))
    require(products > 1e-13, 'Cone did not chemically oxidize fuel')
    report = {'case':'nonuniform-reactive-cone', 'cantera_version':ct.__version__,
              'cells':len(field),'steps':int(history[-1]['steps']), 'final_time':history[-1]['time'],
              'injected_liquid_mass':injected, 'remaining_liquid_mass':remaining,
              'relative_mass_error':mass_error, 'relative_energy_error':energy_error,
              'absolute_momentum_error':momentum_error,'relative_element_error':element_error,
              'absent_element_absolute_kmol':absent_error,
              'temperature_range_K':temperature_range, 'CO_CO2_mass':float(products),
              'restart_fields_parcels_history_exact':True}
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    print('PASS: coupled cone, scheduled injection and exact reactive restart')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',required=True,type=Path)
    parser.add_argument('--yaml',required=True,type=Path)
    parser.add_argument('--case',type=Path,default=Path(__file__).resolve().parents[1]/'cases/spray_3d/methanol_cone.nml')
    parser.add_argument('--output-dir',type=Path)
    args = parser.parse_args()
    if args.output_dir:
        work=args.output_dir.resolve();work.mkdir(parents=True,exist_ok=False)
        check(args.exe.resolve(),args.yaml.resolve(),args.case.resolve(),work)
    else:
        with tempfile.TemporaryDirectory(prefix='pelef-spray-cone-') as directory:
            check(args.exe.resolve(),args.yaml.resolve(),args.case.resolve(),Path(directory))


if __name__ == '__main__':
    main()
