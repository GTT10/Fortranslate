"""Prepare a short, source-bounded decane Spray-Jet application adaptation."""
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
    manifest = json.loads((root/'references/spray_jet_case.json').read_text())
    source = {}
    for name, entry in manifest['files'].items():
        data = (args.reference_source/entry['path']).read_bytes()
        if hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError(f'Pinned jet source mismatch: {name}')
        source[name] = data
    bundle = json.loads(args.bundle.read_text())
    mechanism = json.loads((root/'references/decane_3sp.json').read_text())
    if (not check(args.bundle) or bundle.get('nonreacting') is not True
            or bundle['reactions'] != [] or bundle['species'] != ['NC10H22', 'O2', 'N2']
            or bundle['source']['sha256'] != mechanism['files']['yaml']['sha256']):
        raise ValueError('Expected pinned reaction-free decane bundle')
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
    mw = next(s['molecular_weight'] for s in bundle['thermo'] if s['name'] == 'N2')
    density = 6e6*mw/(8314.46261815324*440)
    fmt = lambda v: format(v, '.17g')
    args.output.mkdir(parents=True, exist_ok=False)
    for name, data in source.items():
        (args.output/('upstream-'+name)).write_bytes(data)
    (args.output/'rate.dat').write_text('# seconds kg/s; all exact upstream knots\n'+
        ''.join(f'{fmt(t)} {fmt(m)}\n' for t, m in points))
    (args.output/'events.dat').write_text('# begin_s end_s port scale\n'+
        ''.join(f'{fmt(a[0])} {fmt(b[0])} 1 1\n' for a, b in intervals))
    (args.output/'level_set.dat').write_text('4 8 4\n'+'1.0\n'*(5*9*5))
    (args.output/'tree.nml').write_text(f'''&eb_tree
 root_shape=4,8,4,child_count=0,
 domain_lower=-.006,0,-.006,domain_upper=.006,.084,.006,
 geometry_type='level_set',level_set_file='level_set.dat',
 time_step={fmt(dt)},final_time={fmt(end_time)},maximum_steps={2*args.events+2},maximum_step_retries=0,
 time_integrator='sdc',sdc_iterations=2,subcycling=.false.,automatic_timestep=.false.,
 chemistry_enabled=.false.,transport_enabled=.true.,simple_transport_enabled=.true.,
 hydro_method='plm_ctu',riemann_solver='pelec',
 boundary_faces='outflow','outflow','no_slip_wall','outflow','outflow','outflow',
 composition_count=1,composition_species='N2',composition_mole_fractions=1,
 regular_density={fmt(density)},cut_density={fmt(density)},initial_pressure=6e6,
 output_file='fortran.csv'
/
''')
    (args.output/'spray.nml').write_text(f'''&eb_spray
 particle_outflow=.true.,particle_options%drag_model='pelephysics',
 particle_options%mixture_time_scheme='euler',
 nc=1,number=4,vapor_names(1)='NC10H22',fractions(1)=1,
 densities(1)=640,heat_capacities(1)=2192.1,
 reference_temperatures(1)=300,reference_pressures(1)=101325,latent_heats(1)=358990,
 mixture_vapor_pressure_models(1)='antoine',
 mixture_antoine_coefficients(:,1)=4.07857,1501.268,-78.67,1e5,
 mixture_antoine_temperature_ranges(:,1)=300,447.27,
 total_mass=1e-10,center=0,1e-8,0,axis=0,1,0,
 speed=1,half_angle={fmt(math.radians(10.5))},diameter=4.5e-6,liquid_temperature=363,
 port_radii(1)=4.5e-5,port_nozzles(1)%profile='solid',
 port_velocity_from_rate(1)=.true.,port_discharge_coefficients(1)=.89,

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
        initial_gas_temperature_K=440, initial_pressure_Pa=6e6,
        event_masses_kg=masses, expected_injected_mass_kg=math.fsum(masses),
        event_speeds_m_s=[m/((b[0]-a[0])*640*math.pi*(4.5e-5)**2*.89) for m, (a, b) in zip(masses, intervals)],
        rate_points=len(points), injection_events=args.events, expected_parcels=4*args.events, final_time_s=end_time,
        adaptations=['4x8x4 all-fluid mesh without AMR', f'first {args.events} active rate intervals',
                     'fuel reference temperature 300 K instead of out-of-range 298.15 K',
                     'nozzle center moved 1e-8 m inside lower-y boundary',
                     'four deterministic parcels per event instead of stochastic sampling',
                     'event-average rate-driven velocity and interval-integrated mass'],
        scope='Application connection; not full Spray-Jet parity or physical validation'), indent=2)+'\n')
    print(args.output.resolve())


if __name__ == '__main__':
    main()
