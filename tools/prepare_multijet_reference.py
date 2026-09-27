"""Prepare a bounded four-port adaptation of the pinned Spray-Multijet case."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from check_mechanism_provenance import check
from pelephysics_transport_import import regenerate_simple_transport


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--reference-source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root/'references/spray_multijet_case.json').read_text())
    source = {}
    for name, entry in manifest['files'].items():
        data = (args.reference_source/entry['path']).read_bytes()
        if hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError(f'Pinned multijet source mismatch: {name}')
        source[name] = data
    bundle = json.loads(args.bundle.read_text())
    mechanism = json.loads((root/'references/decane_3sp.json').read_text())
    if (not check(args.bundle) or bundle.get('nonreacting') is not True
            or bundle['reactions'] != [] or bundle['species'] != ['NC10H22', 'O2', 'N2']
            or bundle['source']['sha256'] != mechanism['files']['yaml']['sha256']):
        raise ValueError('Expected pinned reaction-free decane bundle')
    if bundle.get('simple_transport') != regenerate_simple_transport(bundle, args.bundle.parent):
        raise ValueError('Pinned Simple transport fits required')
    weights = {s['name']: s['molecular_weight'] for s in bundle['thermo']}
    amounts = [.233/weights['O2'], .767/weights['N2']]
    mole = [v/sum(amounts) for v in amounts]
    density = 2.93e6/(8314.46261815324*900*sum(amounts))
    dt = 1e-7
    centers = [(x, 1e-8, z) for x in (.003, .009) for z in (.003, .009)]
    fmt = lambda v: format(v, '.17g')
    args.output.mkdir(parents=True, exist_ok=False)
    for name, data in source.items():
        (args.output/('upstream-'+name)).write_bytes(data)
    (args.output/'level_set.dat').write_text('4 8 4\n'+'1.0\n'*(5*9*5))
    (args.output/'events.dat').write_text('# begin_s end_s port scale\n'+''.join(
        f'{fmt(step*dt)} {fmt((step+1)*dt)} {port} 1\n'
        for step in range(2) for port in range(1, 5)))
    (args.output/'rate.dat').write_text('0 .0023\n1.54e-3 .0023\n')
    (args.output/'tree.nml').write_text(f'''&eb_tree
 root_shape=4,8,4,child_count=0,
 domain_lower=0,0,0,domain_upper=.012,.012,.012,
 geometry_type='level_set',level_set_file='level_set.dat',
 time_step={fmt(dt)},final_time={fmt(2*dt)},maximum_steps=4,maximum_step_retries=0,
 time_integrator='sdc',sdc_iterations=2,subcycling=.false.,automatic_timestep=.false.,
 chemistry_enabled=.false.,transport_enabled=.true.,simple_transport_enabled=.true.,
 hydro_method='plm_ctu',riemann_solver='pelec',
 boundary_faces='periodic','periodic','no_slip_wall','outflow','periodic','periodic',
 composition_count=2,composition_species='O2','N2',
 composition_mole_fractions={fmt(mole[0])},{fmt(mole[1])},
 regular_density={fmt(density)},cut_density={fmt(density)},initial_pressure=2.93e6,
 output_file='fortran.csv'
/
''')
    ports = ''.join(f' port_centers(:,{i})='+','.join(map(fmt, center))+',\n'
                    for i, center in enumerate(centers, 1))
    (args.output/'spray.nml').write_text(f'''&eb_spray
 particle_outflow=.true.,particle_options%drag_model='pelephysics',
 particle_options%mixture_time_scheme='euler',
 nc=1,number=4,vapor_names(1)='NC10H22',fractions(1)=1,
 densities(1)=669,heat_capacities(1)=2593.9,
 reference_temperatures(1)=300,reference_pressures(1)=101325,latent_heats(1)=302940,
 mixture_vapor_pressure_models(1)='antoine',
 mixture_antoine_coefficients(:,1)=4.10549,1625.928,-92.839,1e5,
 mixture_antoine_temperature_ranges(:,1)=300,489,
 total_mass=1e-10,axis=0,1,0,speed=600,
 half_angle={fmt(math.radians(10.5))},diameter=4.5e-6,liquid_temperature=363,
 port_count=4,port_radii(1:4)=4*4.5e-5,
 port_nozzles(1:4)%profile='solid','solid','solid','solid',
{ports}/
&eb_spray_rate_files
 profile_files(1:4)='rate.dat','rate.dat','rate.dat','rate.dat'
/
&eb_spray_events_file
 event_file='events.dat'
/
''')
    (args.output/'preparation.json').write_text(json.dumps(dict(
        reference=manifest, bundle_sha256=hashlib.sha256(args.bundle.read_bytes()).hexdigest(),
        initial_gas_temperature_K=900, initial_pressure_Pa=2.93e6,
        initial_mass_fractions=dict(O2=.233, N2=.767), port_centers_m=centers,
        mass_flow_per_port_kg_s=.0023, event_masses_kg=[.0023*dt]*8,
        expected_injected_mass_kg=8*.0023*dt, injection_events=8,
        expected_parcels=32, final_time_s=2*dt,
        adaptations=['4x8x4 all-fluid mesh without AMR', 'two 100 ns intervals per port',
                     'fuel reference temperature 300 K instead of out-of-range 298.15 K',
                     'nozzle centers moved 1e-8 m inside lower-y boundary',
                     'four deterministic parcels per event instead of stochastic sampling',
                     'Antoine evaluation restricted to 300-489 K; not a fitted validity claim'],
        scope='Four-port application connection; not full Spray-Multijet parity or physical validation'),
        indent=2)+'\n')
    print(args.output.resolve())


if __name__ == '__main__':
    main()
