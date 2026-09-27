"""Prepare pinned Spray-EB geometry and five physical initial droplets in SI units."""
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
    parser.add_argument('--cells', type=int, default=64)
    parser.add_argument('--final-time', type=float, default=0.000273365690540257)
    args = parser.parse_args()
    if not 4 <= args.cells <= 128 or not math.isfinite(args.final_time) or args.final_time <= 0:
        raise ValueError('Require 4..128 cells and a positive finite final time')
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root/'references/spray_eb_case.json').read_text())
    source = {}
    for name, entry in manifest['files'].items():
        data = (args.reference_source/entry['path']).read_bytes()
        if hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError(f'Pinned Spray-EB source mismatch: {name}')
        source[name] = data
    bundle = json.loads(args.bundle.read_text())
    mechanism = json.loads((root/'references/heptane_3sp.json').read_text())
    if (not check(args.bundle) or bundle.get('nonreacting') is not True
            or bundle['reactions'] != [] or bundle['species'] != ['NC7H16', 'O2', 'N2']
            or bundle['source']['sha256'] != mechanism['files']['yaml']['sha256']):
        raise ValueError('Expected pinned reaction-free heptane bundle')
    if bundle.get('simple_transport') != regenerate_simple_transport(bundle, args.bundle.parent):
        raise ValueError('Pinned Simple transport fits required')
    lines = source['initspraydata_3d.dat'].decode().splitlines()
    if int(lines[0]) != 5 or len(lines) != 6:
        raise ValueError('Expected five upstream initial particles')
    particles = []
    for line in lines[1:]:
        row = list(map(float, line.split()))
        # xyz, velocity xyz, T, diameter, Yfuel, number, N0, BM1, BM2, film height.
        if len(row) != 14 or row[8:] != [1, 1, 1, 0, 0, 0]:
            raise ValueError('Unsupported upstream initial particle state')
        diameter = row[7]*.01
        velocity = [v*.01 for v in row[3:6]]
        speed = math.sqrt(sum(v*v for v in velocity))
        particles.append(dict(position_m=[v*.01 for v in row[:3]], velocity_m_s=velocity,
                              axis=[v/speed for v in velocity], speed_m_s=speed,
                              temperature_K=row[6], diameter_m=diameter,
                              mass_kg=693*math.pi/6*diameter**3))
    weights = {s['name']: s['molecular_weight'] for s in bundle['thermo']}
    amounts = [.233/weights['O2'], .767/weights['N2']]
    mole = [v/sum(amounts) for v in amounts]
    density = 1.01325e6/(8314.46261815324*1500*sum(amounts))
    n = args.cells
    point = [.08830222, .01786062, 0]
    normal = [.7660444, -.6427876, 0]
    fmt = lambda v: format(v, '.17g')
    args.output.mkdir(parents=True, exist_ok=False)
    for name, data in source.items():
        (args.output/('upstream-'+name)).write_bytes(data)
    with (args.output/'level_set.dat').open('w') as out:
        out.write(f'{n} {n} {n}\n')
        for k in range(n+1):
            for j in range(n+1):
                for i in range(n+1):
                    # AMReX plane IF is negative in fluid; PeleF is positive.
                    phi = -sum((x-point[d])*normal[d] for d, x in enumerate([.1*i/n, .1*j/n, .1*k/n]))
                    out.write(fmt(phi)+'\n')
    (args.output/'tree.nml').write_text(f'''&eb_tree
 root_shape={n},{n},{n},child_count=0,
 domain_lower=0,0,0,domain_upper=.1,.1,.1,
 geometry_type='level_set',level_set_file='level_set.dat',
 time_step=1e-7,final_time={fmt(args.final_time)},maximum_steps=100000,maximum_step_retries=12,
 time_integrator='sdc',sdc_iterations=2,automatic_timestep=.true.,hydro_cfl=.3,
 chemistry_enabled=.false.,transport_enabled=.true.,simple_transport_enabled=.true.,
 hydro_method='plm_ctu',riemann_solver='pelec',redistribution='state_redist',
 boundary_faces=4*'no_slip_wall',2*'periodic',
 embedded_wall_temperature=1500,embedded_noslip=.true.,
 composition_count=2,composition_species='O2','N2',
 composition_mole_fractions={fmt(mole[0])},{fmt(mole[1])},
 regular_density={fmt(density)},cut_density={fmt(density)},initial_pressure=1.01325e6,
 output_file='fortran.csv'
/
''')
    ports = ''
    for i, p in enumerate(particles, 1):
        ports += f' port_centers(:,{i})='+','.join(map(fmt, p['position_m']))+',\n'
        ports += f' port_axes(:,{i})='+','.join(map(fmt, p['axis']))+',\n'
        for field, key in [('masses', 'mass_kg'), ('speeds', 'speed_m_s'), ('diameters', 'diameter_m')]:
            ports += f' port_{field}({i})={fmt(p[key])},\n'
    (args.output/'spray.nml').write_text(f'''&eb_spray
 particle_outflow=.true.,particle_options%drag_model='pelephysics',
 particle_options%mixture_time_scheme='euler',
 nc=1,number=1,vapor_names(1)='NC7H16',fractions(1)=1,
 densities(1)=693,heat_capacities(1)=2248.3,
 reference_temperatures(1)=300,reference_pressures(1)=101325,latent_heats(1)=363000,
 mixture_vapor_pressure_models(1)='antoine',
 mixture_antoine_coefficients(:,1)=4.02832,1268.636,-56.199,1e5,
 mixture_antoine_temperature_ranges(:,1)=300,371.6,
 half_angle=0,liquid_temperature=300,
 port_count=5,injection_events=5,event_ports(1:5)=1,2,3,4,5,
 injection_times(1:5)=5*0,
{ports}/
''')
    adaptations = ['Initial particles represented by five instantaneous events at t=0',
                   'CTU/SDC instead of upstream method-of-lines integration',
                   'Single root patch instead of upstream max_grid_size=32 partition',
                   'Time step capped at 100 ns with adaptive CFL and retry controls',
                   'Antoine evaluation restricted to 300-371.6 K; not a fitted validity claim',
                   'Particle wall interactions require separate reference comparison']
    if n != 64:
        adaptations.append(f'{n} cubed grid instead of upstream 64 cubed')
    if args.final_time != 0.000273365690540257:
        adaptations.append('Changed final time for connection check')
    (args.output/'preparation.json').write_text(json.dumps(dict(
        reference=manifest, bundle_sha256=hashlib.sha256(args.bundle.read_bytes()).hexdigest(),
        cells=n, final_time_s=args.final_time, initial_gas_temperature_K=1500,
        initial_pressure_Pa=1.01325e6, initial_mass_fractions=dict(O2=.233, N2=.767),
        plane_point_m=point, upstream_plane_normal=normal, particles=particles,
        expected_injected_mass_kg=sum(p['mass_kg'] for p in particles),
        adaptations=adaptations, scope='Source-preserving case connection; full solver parity remains unproven'),
        indent=2)+'\n')
    print(args.output.resolve())


if __name__ == '__main__':
    main()
