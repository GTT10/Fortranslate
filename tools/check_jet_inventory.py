"""Check completion and fuel accounting for the prepared single-fuel decane jet."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re


def rows(path):
    with path.open() as stream:
        result = list(csv.DictReader(stream))
    if any(not math.isfinite(float(v)) for row in result for v in row.values()):
        raise ValueError(f'Nonfinite output: {path}')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    folder = args.run_dir
    metadata = json.loads((folder/'preparation.json').read_text())
    log = (folder/'fortran.log').read_text()
    complete = re.findall(r'Completed steps:\s*(\d+)\s*time:\s*([\d.Ee+-]+)', log)
    events = re.findall(r'Processed injection events:\s*(\d+)', log)
    if not complete or not events:
        raise ValueError('No completed application output')
    expected_events = metadata.get('injection_events', len(metadata['event_masses_kg']))
    if int(events[-1]) != expected_events or not math.isclose(
            float(complete[-1][1]), metadata['final_time_s'], rel_tol=1e-13, abs_tol=0):
        raise ValueError('Application did not reach the prepared schedule/time')
    files = [folder/'fortran.csv', folder/'fortran.csv.budget.csv']
    particle_files = sorted(folder.glob('fortran.csv.particles.rank*.csv'))
    film_files = sorted(folder.glob('fortran.csv.films.rank*.csv'))
    if not particle_files or not film_files:
        raise ValueError('Missing particle or film output')
    files += particle_files+film_files
    gas = rows(files[0])
    budget = {int(v['component']): v for v in rows(files[1])}
    particles = [v for path in particle_files for v in rows(path)]
    if any(rows(path) for path in film_files):
        raise ValueError('Unexpected film in the free-jet fixture')
    if len(gas) != 128 or any(float(v['fluid_volume']) <= 0 for v in gas):
        raise ValueError('Expected 128 all-fluid cells')
    keys = [(v['node'], v['level'], v['i'], v['j'], v['k']) for v in gas]
    ids = [int(v['id']) for v in particles]
    if len(set(keys)) != len(gas) or len(set(ids)) != len(ids):
        raise ValueError('Duplicate cells or parcel IDs')
    if len(particles) > metadata['expected_parcels'] or any(i <= 0 for i in ids):
        raise ValueError('Unexpected parcel population')
    for p in particles:
        if float(p['mass']) <= 0 or float(p['multiplicity']) <= 0 or not math.isclose(
                float(p['component_mass_1']), float(p['mass']), rel_tol=1e-12):
            raise ValueError('Invalid single-fuel parcel')
    expected = metadata['expected_injected_mass_kg']
    if expected <= 0 or not math.isfinite(expected):
        raise ValueError('Invalid prepared fuel mass')
    vapor = math.fsum(float(v['fluid_volume'])*float(v['rhoY_NC10H22']) for v in gas)
    liquid = math.fsum(float(v['mass'])*float(v['multiplicity']) for v in particles)
    fuel = {k: float(v) for k, v in budget[6].items()}
    if any(fuel[k] != 0 for k in ('initial', 'chemistry', 'force')):
        raise ValueError('Unexpected initial/reaction/force fuel contribution')
    # The application ledger records signed additions; outflow is negative.
    reconstructed = vapor+liquid-fuel['boundary']-fuel['wall']
    error = abs(reconstructed-expected)/expected
    inventory_error = abs(vapor+liquid-fuel['final'])/expected
    injected_error = abs(fuel['injected']-expected)/expected
    result = dict(completed_steps=int(complete[-1][0]), completed_time_s=float(complete[-1][1]),
                  processed_events=int(events[-1]), retained_parcels=len(particles),
                  expected_injected_mass_kg=expected, gas_fuel_mass_kg=vapor, liquid_fuel_mass_kg=liquid,
                  signed_boundary_fuel_kg=fuel['boundary'], signed_wall_fuel_kg=fuel['wall'],
                  fuel_balance_relative_error=error, csv_inventory_relative_error=inventory_error,
                  injected_mass_relative_error=injected_error, tolerance=1e-8,
                  energy_budget_residual_J=float(budget[5]['residual']),
                  sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in files+[folder/'preparation.json', folder/'fortran.log']})
    result['passed'] = max(error, inventory_error, injected_error) <= 1e-8
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    if not result['passed']:
        raise ValueError('Jet fuel accounting failed; see report')
    print(f"PASS: {result['processed_events']} events, {len(particles)} retained parcels, fuel error {error:.3e}")


if __name__ == '__main__':
    main()
