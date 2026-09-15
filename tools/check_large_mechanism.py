#!/usr/bin/env python3
"""Exercise a pinned GRI53 or Lu n-dodecane bundle through selected 0D and spray executables."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    import cantera as ct
    parser = argparse.ArgumentParser()
    parser.add_argument('--mechanism', choices=['gri53', 'dodecane'], default='gri53')
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    dodecane = args.mechanism == 'dodecane'
    work = Path(tempfile.mkdtemp(prefix=args.mechanism+'-', dir=args.work_dir.resolve()))
    def run(executable, name, text):
        path = work / (name + '.nml')
        path.write_text(text.rstrip()+"\n")
        result = subprocess.run([str(args.build_dir.resolve() / executable), str(path)], cwd=work,
                                capture_output=True, text=True, timeout=180)
        (work / (name + '.log')).write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
    reactor = '''&selected_reactor
 integrator='cvode', initial_temperature=1500, initial_pressure=101325,
 composition_count=3,composition_species='CH4','O2','N2',composition_mole_fractions=1,2,7.52,
 final_time=1e-4,output_interval=2.5e-5,initial_time_step=1e-12,
 minimum_time_step=1e-16,maximum_time_step=1e-6,relative_tolerance=1e-9,absolute_tolerance=1e-15,
 maximum_steps=100000,output_file='reactor.csv'
/'''
    if dodecane:
        reactor = reactor.replace('1500', '1100').replace('101325', '4000000').replace("'CH4'", "'NC12H26'")
        reactor = reactor.replace('1,2,7.52', '1,37,139.12').replace('final_time=1e-4', 'final_time=6e-4')
        reactor = reactor.replace('output_interval=2.5e-5', 'output_interval=1e-5')
    run('pelef0d_selected', 'reactor', reactor)
    gas = ct.Solution(str(args.source.resolve()), 'gas' if dodecane else 'gri30')
    assert gas.n_species == 53 and gas.n_reactions == (268 if dodecane else 325)
    gas.TPX = ((1100, 4e6, 'NC12H26:1,O2:37,N2:139.12') if dodecane else
               (1500, 101325, 'CH4:1,O2:2,N2:7.52'))
    cv = ct.IdealGasReactor(gas, energy='on', clone=True)
    network = ct.ReactorNet([cv]); network.rtol=1e-11;network.atol=1e-18
    max_t=max_y=0.0
    with (work / 'reactor.csv').open() as stream:
        records = list(csv.DictReader(stream))
    for row in records:
        time=float(row['time'])
        if time: network.advance(time)
        max_t=max(max_t,abs(float(row['temperature'])-cv.T))
        max_y=max(max_y,max(abs(float(row['Y_'+name])-cv.phase.Y[k]) for k,name in enumerate(gas.species_names)))
    assert max_t < 1e-3 and max_y < 1e-7, (max_t,max_y)
    if dodecane:
        assert float(records[-1]['temperature']) > 2000, 'Fixture must ignite'
    fuel_column = 'Y_NC12H26' if dodecane else 'Y_CH4'
    assert abs(float(records[-1][fuel_column])-float(records[0][fuel_column])) > 1e-7
    base = '''&spray
 cells=4,4,4,domain_length=.004,.004,.004,
 gas_temperature=1500,gas_pressure=101325,mole_names='CH4','O2','N2',mole_amounts=1,2,7.52,
 chemistry_enabled=.true.,chemistry_method='cvode',transport_enabled=.true.,
 vapor_name='H2O',liquid_density=1000,liquid_cp=4184,
 reference_temperature=373.15,reference_pressure=101325,latent_heat=2257000,
 drop_temperature=330,drop_diameter=20e-6,initial_liquid_mass=1e-10,parcel_count=2,
 injection_center=.002,.002,.002,injection_axis=0,0,1,injection_speed=10,
 final_time=1e-7,maximum_dt=5e-8,
 output_prefix='TAG'
/'''
    if dodecane:
        base = base.replace("'CH4'", "'NC12H26'").replace('1,2,7.52', '1,37,139.12')
    for name,extra in [('reference',''),('stop',"stop_after_steps=1,checkpoint_file='spray.chk',"),
                       ('resumed',"restart_file='spray.chk',")]:
        run('pelef_spray_3d_selected',name,base.replace('TAG',name).replace('&spray','&spray\n'+extra))
    for suffix in ['gas.csv','parcels.csv']:
        assert (work / ('reference.'+suffix)).read_bytes() == (work / ('resumed.'+suffix)).read_bytes()
    def last(name):
        with (work / (name+'.history.csv')).open() as stream:
            return list(csv.DictReader(stream))[-1]
    row=last('reference');assert row==last('resumed')
    assert abs(float(row['energy_error'])) < 1e-9*abs(float(row['total_energy']))
    result={'mechanism': args.mechanism, 'final_temperature_K': float(records[-1]['temperature']),
            'species':gas.n_species,'reactions':gas.n_reactions,'maximum_temperature_error_K':max_t,
            'maximum_mass_fraction_error':max_y,'spray_exact_restart':True,'artifacts':str(work)}
    (work/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
