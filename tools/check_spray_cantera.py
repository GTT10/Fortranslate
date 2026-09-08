#!/usr/bin/env python3
"""Independent Cantera chemistry checks and a closed methanol spray burn gate.

Chemistry-only cases compare complete temperature histories and final fields
with Cantera. The spray case checks conservation, elements and actual fuel
conversion; it is an integration test, not experimental spray validation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import cantera as ct
import numpy as np


def rows(path):
    with path.open() as stream:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def run(args, work):
    require(ct.__version__ == '3.2.0', 'Reference requires Cantera 3.2.0')
    spray = args.case == 'methanol-spray'
    methanol = args.case.startswith('methanol')
    names = ['O2','N2',''] if spray else ['CH3OH','O2','N2'] if methanol else ['H2','O2','N2']
    amounts = [1.,3.76,0.] if spray else [1.,1.5,5.64] if methanol else [2.,1.,3.76]
    t0 = 2200. if spray else 1800. if methanol else 1400.
    pressure = 5e5 if methanol else 101325.
    end = 5e-5 if spray else 3e-5 if methanol else 2e-5
    dt = args.maximum_dt or 2e-7
    liquid_mass = 1.6e-9 if spray else 0.
    prefix = work/'run'
    material = ''
    if methanol:
        material = """vapor_name='CH3OH',liquid_density=792,liquid_cp=2530,
 reference_temperature=337.85,reference_pressure=101325,latent_heat=1100000,"""
    inp = work/'input.nml'
    inp.write_text(f'''&spray
 cells=2,2,2,domain_length=.004,.004,.004,hydro_reconstruction='pcm',
 gas_temperature={t0},gas_pressure={pressure},
 mole_names={','.join(repr(n) for n in names)},mole_amounts={','.join(str(a) for a in amounts)},
 initial_liquid_mass={liquid_mass},parcel_count=1,injection_center=.002,.002,.002,
 injection_speed=0,cone_angle=0,drop_diameter=2e-6,drop_temperature=300,
 final_time={end},maximum_dt={dt},chemistry_method='{args.backend}',
 rtol=1e-7,atol=1e-13,transport_enabled=.false.,smagorinsky_constant=0,
 {material}
 output_prefix='{prefix}',checkpoint_file='{work/'final.chk'}'
/\n''')
    result = subprocess.run([str(args.exe.resolve()), str(inp)], cwd=work, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
    (work/'stdout.log').write_text(result.stdout)
    require(result.returncode == 0, result.stdout)
    history = rows(work/'run.history.csv')
    field = rows(work/'run.gas.csv')
    require(len(field)==8, 'Expected eight periodic cells')
    require(history[-1]['time'] == end, 'Requested end time not reached')
    require(all(np.isfinite(v) for row in history+field for v in row.values()), 'Nonfinite output')
    mass0 = history[0]['system_mass']; energy0=history[0]['total_energy']
    mass_error=max(abs(r['mass_error']) for r in history)/mass0
    energy_error=max(abs(r['energy_error']) for r in history)/abs(energy0)
    require(mass_error < 2e-11 and energy_error < 2e-11, 'Global mass/energy conservation failed')
    phase = ct.Solution(str(args.yaml.resolve()), args.phase)
    phase.TPX=t0,pressure,{n:a for n,a in zip(names,amounts) if n}
    rho0=phase.density; volume=.004**3; gas_mass0=rho0*volume
    matrix=np.array([[phase.n_atoms(k,el)/phase.molecular_weights[k] for k in range(phase.n_species)]
                     for el in phase.element_names])
    initial_moles=matrix @ (phase.Y*gas_mass0)
    if spray:
        initial_moles+=matrix[:,phase.species_index('CH3OH')]*liquid_mass
    mean_species=np.mean([[r['rhoY_'+n] for n in phase.species_names] for r in field],axis=0)*volume
    parcel_data=rows(work/'run.parcels.csv')
    remaining=sum(p['mass']*p['multiplicity'] for p in parcel_data)
    final_moles=matrix @ mean_species
    if spray:
        final_moles+=matrix[:,phase.species_index('CH3OH')]*remaining
    element_error=float(np.max(np.abs(final_moles-initial_moles)/np.maximum(initial_moles,1e-25)))
    require(element_error<2e-8, f'Element inventory mismatch: {element_error}')
    report={'case':args.case,'backend':args.backend,'cantera_version':ct.__version__,
            'source_sha256':hashlib.sha256(args.yaml.read_bytes()).hexdigest(),
            'steps':int(history[-1]['steps']),'final_time':end,
            'relative_mass_error':mass_error,'relative_energy_error':energy_error,
            'relative_element_error':element_error,'remaining_liquid_mass':remaining,
            'initial_temperature':t0,'final_temperature':field[0]['temperature']}
    if spray:
        carbon_to_oxides=sum(mean_species[phase.species_index(n)]/phase.molecular_weights[phase.species_index(n)]
                             for n in ['CO','CO2'])/(liquid_mass/phase.molecular_weights[phase.species_index('CH3OH')])
        report['carbon_in_CO_and_CO2_fraction']=float(carbon_to_oxides)
        require(remaining<.01*liquid_mass, 'Spray did not substantially evaporate')
        require(carbon_to_oxides>.5, 'No substantial oxidation of injected liquid fuel')
        require(field[0]['temperature']>t0+50, 'No net heat release after evaporation')
    else:
        reactor=ct.IdealGasReactor(phase,clone=True)
        network=ct.ReactorNet([reactor]); network.rtol=1e-11; network.atol=1e-18
        errors=[]
        for row in history:
            if row['time']>network.time: network.advance(row['time'])
            errors.append(max(abs(row['tmin']-reactor.T),abs(row['tmax']-reactor.T)))
        ys=np.array([[r['rhoY_'+n]/r['rho'] for n in phase.species_names] for r in field])
        y_error=float(np.max(np.abs(ys-reactor.phase.Y)))
        report.update(maximum_temperature_error_K=max(errors),maximum_mass_fraction_error=y_error,
                      cantera_final_temperature=reactor.T)
        require(max(errors)<.08, f'Cantera temperature mismatch: {max(errors)} K')
        require(y_error<3e-6, f'Cantera species mismatch: {y_error}')
        require(reactor.T>t0+100, 'Reference did not traverse ignition')
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    print('PASS: '+args.case)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',required=True,type=Path)
    parser.add_argument('--yaml',required=True,type=Path)
    parser.add_argument('--phase',required=True)
    parser.add_argument('--case',choices=['h2-chemistry','methanol-chemistry','methanol-spray'],required=True)
    parser.add_argument('--backend',choices=['implicit','cvode'],default='cvode')
    parser.add_argument('--maximum-dt',type=float)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    if args.maximum_dt is not None and (not np.isfinite(args.maximum_dt) or args.maximum_dt<=0):
        parser.error('--maximum-dt must be positive and finite')
    if args.output_dir:
        work=args.output_dir.resolve();work.mkdir(parents=True,exist_ok=False);run(args,work)
    else:
        with tempfile.TemporaryDirectory(prefix='pelef-spray-reference-') as directory:
            run(args,Path(directory))


if __name__=='__main__':
    main()
