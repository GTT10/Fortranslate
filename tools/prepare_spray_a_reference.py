"""Prepare a source-bounded n-dodecane Spray-A/KHRT application adaptation."""
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
    parser.add_argument('--events', type=int, default=2, help='Number of active ROI intervals to prepare')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root/'references/spray_a_case.json').read_text())
    source = {}
    for name, entry in manifest['files'].items():
        data = (args.reference_source/entry['path']).read_bytes()
        if hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError(f'Pinned Spray-A source mismatch: {name}')
        source[name] = data
    bundle = json.loads(args.bundle.read_text())
    mechanism = json.loads((root/'references/dodecane_lu.json').read_text())
    if (not check(args.bundle) or len(bundle['species']) != 53
            or len(bundle['reactions']) != 268 or bundle['species'][0] != 'NC12H26'
            or bundle['source']['sha256'] != mechanism['files']['yaml']['sha256']):
        raise ValueError('Expected pinned Lu n-dodecane bundle')
    if bundle.get('simple_transport') != regenerate_simple_transport(bundle, args.bundle.parent):
        raise ValueError('Pinned Simple transport fits required')
    lines = source['roi.dat'].decode('ascii').splitlines()
    if lines[0] != 'Time [ms];Mass flow rate [g/s]':
        raise ValueError('Unexpected rate units')
    points = [tuple(float(v)*1e-3 for v in line.split(';')) for line in lines[1:]]
    if (any(len(row) != 2 or not all(math.isfinite(v) for v in row)
            or min(row) < 0 for row in points)
            or any(b[0] <= a[0] for a, b in zip(points, points[1:]))):
        raise ValueError('Invalid reference rate table')
    # First positive knot is the pinned SprayJet::readROI start time.
    if points[0] != (0., 0.) or points[1][1] <= 0:
        raise ValueError('Unexpected injection start')
    intervals = [(a, b) for a, b in zip(points[1:], points[2:]) if a[1] > 0 or b[1] > 0]
    if not 1 <= args.events <= len(intervals):
        raise ValueError(f'Event count must be between 1 and {len(intervals)}')
    intervals = intervals[:args.events]
    dt = points[1][0]
    end_time = intervals[-1][1][0]
    masses = [.5*(b[0]-a[0])*(a[1]+b[1]) for a, b in intervals]
    composition = dict(N2=.8971, H2O=.0377, CO2=.0652)
    weights = {s['name']: s['molecular_weight'] for s in bundle['thermo']}
    mw = sum(x*weights[name] for name,x in composition.items())
    density = 6e6*mw/(8314.46261815324*900)
    rho_coefficients = [1779.9145710513213, -6.848424883513244, .015016559178741312, -1.2149370117669905e-5]
    liquid_density = sum(a*363**i for i,a in enumerate(rho_coefficients))
    fmt = lambda v: format(v, '.17g')
    args.output.mkdir(parents=True, exist_ok=False)
    for name, data in source.items():
        (args.output/('upstream-'+name)).write_bytes(data)
    (args.output/'rate.dat').write_text('# seconds kg/s; all exact upstream knots\n'+
        ''.join(f'{fmt(t)} {fmt(m)}\n' for t, m in points))
    (args.output/'events.dat').write_text('# begin_s end_s port scale\n'+
        ''.join(f'{fmt(a[0])} {fmt(b[0])} 1 1\n' for a, b in intervals))
    (args.output/'level_set.dat').write_text('4 4 8\n'+'1.0\n'*(5*5*9))
    (args.output/'tree.nml').write_text(f'''&eb_tree
 root_shape=4,4,8,child_count=0,
 domain_lower=-.006,-.006,0,domain_upper=.006,.006,.06,
 geometry_type='level_set',level_set_file='level_set.dat',
 time_step={fmt(dt)},final_time={fmt(end_time)},maximum_steps={math.ceil(end_time/dt)+args.events+2},maximum_step_retries=0,
 time_integrator='sdc',sdc_iterations=2,subcycling=.false.,automatic_timestep=.false.,
 chemistry_enabled=.false.,transport_enabled=.true.,simple_transport_enabled=.true.,
 hydro_method='plm_ctu',riemann_solver='pelec',
 boundary_faces='outflow','outflow','outflow','outflow','no_slip_wall','outflow',
 composition_count=3,composition_species='N2','H2O','CO2',composition_mole_fractions=.8971,.0377,.0652,
 regular_density={fmt(density)},cut_density={fmt(density)},initial_pressure=6e6,
 output_file='fortran.csv'
/
''')
    (args.output/'spray.nml').write_text(f'''&eb_spray
 particle_outflow=.true.,particle_options%drag_model='pelephysics',
 particle_options%mixture_time_scheme='euler',
 nc=1,number=4,vapor_names(1)='NC12H26',fractions(1)=1,
 densities(1)={fmt(liquid_density)},heat_capacities(1)=2217,
 mixture_density_models(1)='polynomial',
 mixture_density_coefficients(:,1)={','.join(map(fmt,rho_coefficients))},
 mixture_density_temperature_ranges(:,1)=300,658,
 reference_temperatures(1)=300,reference_pressures(1)=101325,latent_heats(1)=359411,
 mixture_boiling_models(1)='pelemp-cap',
 mixture_critical_temperatures(1)=658,mixture_normal_boiling_temperatures(1)=489,
 mixture_molecular_weights(1)={fmt(weights['NC12H26'])},
 mixture_vapor_pressure_models(1)='antoine',
 mixture_antoine_coefficients(:,1)=4.265416493004907,1741.1453093807672,-80.90351868881419,1e5,
 mixture_antoine_temperature_ranges(:,1)=300,658,
 total_mass=1e-10,center=0,0,1e-8,axis=0,0,1,
 speed=1,half_angle={fmt(math.radians(10))},diameter=89e-6,liquid_temperature=363,
 port_radii(1)=4.5e-5,port_nozzles(1)%profile='solid',
 port_velocity_from_rate(1)=.true.,port_discharge_coefficients(1)=.89,
 particle_options%breakup_model='khrt',particle_options%surface_tension=.015,
 particle_options%khrt_b0=.61,particle_options%khrt_b1=7.3,particle_options%khrt_c3=5.3,
 particle_options%khrt_children=.true.,particle_options%khrt_parcel_factor=.5,
 particle_options%khrt_reference_diameter_cubed={fmt((89e-6)**3)},
 particle_options%liquid_viscosity_model='inverse-temperature',
 particle_options%liquid_viscosity_coefficients=-.0028611188933370675,3.940765272136905,-1812.5486152106896,302158.7187559424,
 particle_options%liquid_viscosity_temperature_range=300,658,

/
&eb_spray_rate_files
 profile_files(1)='rate.dat'
/
&eb_spray_events_file
 event_file='events.dat'
/
''')
    (args.output/'preparation.json').write_text(json.dumps(dict(
        reference=manifest, bundle_sha256=hashlib.sha256(args.bundle.read_bytes()).hexdigest(),
        initial_gas_temperature_K=900, initial_mole_fractions=composition, liquid_density_at_injection_kg_m3=liquid_density, breakup_model="KHRT", chemistry_enabled=False, initial_pressure_Pa=6e6,
        event_masses_kg=masses, expected_injected_mass_kg=math.fsum(masses),
        event_speeds_m_s=[m/((b[0]-a[0])*liquid_density*math.pi*(4.5e-5)**2*.89) for m, (a, b) in zip(masses, intervals)],
        rate_points=len(points), injection_events=args.events, expected_parcels=4*args.events, final_time_s=end_time,
        adaptations=['4x4x8 all-fluid mesh without AMR', f'first {args.events} active rate intervals',
                     'fuel reference temperature 300 K instead of out-of-range 298.15 K',
                     'nozzle center moved 1e-8 m inside lower-z boundary',
                     'four deterministic parcels per event instead of stochastic sampling',
                     'event-average rate-driven velocity and interval-integrated mass',
                     'liquid property evaluation restricted to 300-658 K; not a fit-validity claim'],
        scope='Application connection; not full Spray-A parity or physical validation'), indent=2)+'\n')
    print(args.output.resolve())


if __name__ == '__main__':
    main()
