#!/usr/bin/env python3
"""Generate an explicit LiDryer variant with the pinned native activation basis.

The canonical Cantera/SI mechanism is untouched. Only activation temperatures
are changed, using the pinned header's own calorie gas constant. This is a
diagnostic compatibility variant, not a new physical mechanism or qualification.
Run with the same Cantera environment as the canonical bundle importer.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def equation_key(equation):
    """Canonical stoichiometry for the pinned LiDryer equation syntax."""
    sides = re.sub(r'\(\s*\+\s*M\s*\)', '+ M', equation).split('<=>')
    if len(sides) != 2:
        raise ValueError('Unsupported native reaction equation')
    result = []
    for side in sides:
        counts = {}
        for term in side.split('+'):
            match = re.fullmatch(r'\s*(?:(\d+)\s+)?([A-Za-z][A-Za-z0-9]*)\s*', term)
            if match is None:
                raise ValueError('Unsupported native reaction term')
            count, name = match.groups()
            counts[name] = counts.get(name, 0) + int(count or 1)
        result.append(sorted(counts.items()))
    return result


def activation_mapping(source, header, rc):
    """Match each pinned reaction and pressure branch, including zero Ea."""
    number = r'[-+0-9.eE]+'
    reactions = re.findall(r'^- equation:.*?(?=^- equation:|\Z)', source,
                           flags=re.MULTILINE | re.DOTALL)
    section = header[header.index('comp_qfqr('):header.index('productionRate(')]
    section = section[section.index('// Evaluate the kfs'):]
    blocks = re.findall(r'// reaction (\d+):([^\n]*)\n(.*?)(?=// reaction |\Z)',
                        section, flags=re.DOTALL)
    if len(reactions) != 21 or [int(b[0]) for b in blocks] != list(range(21)):
        raise ValueError('Unexpected LiDryer reaction mapping')
    mapping = []
    for index, (source_reaction, (_, equation, block)) in enumerate(zip(reactions, blocks)):
        source_equation = source_reaction.splitlines()[0].split(':', 1)[1].split('#', 1)[0].strip()
        if equation_key(source_equation) != equation_key(equation):
            raise ValueError(f'Native reaction equation mismatch at reaction {index + 1}')
        rates = re.findall(
            r'^  (rate-constant|high-P-rate-constant|low-P-rate-constant):'
            r'[^\n]*\bEa:\s*(' + number + r')', source_reaction, flags=re.MULTILINE)
        names = [name for name, _ in rates]
        if sorted(names) not in (['rate-constant'], ['high-P-rate-constant', 'low-P-rate-constant']):
            raise ValueError(f'Unexpected source rate branches at reaction {index + 1}')
        if bool(re.search(r'\bredP\s*=', block)) != ('low-P-rate-constant' in names):
            raise ValueError(f'Native pressure branches differ at reaction {index + 1}')
        for branch, energy in rates:
            variable = 'redP' if branch == 'low-P-rate-constant' else 'k_f'
            assignments = re.findall(r'\b' + variable + r'\s*=\s*([^;]+);', block)
            if len(assignments) != 1:
                raise ValueError(f'Invalid native {variable} at reaction {index + 1}')
            expression = assignments[0]
            values = re.findall(r'-\s*\((' + number + r')\)\s*\*\s*invT', expression)
            if len(values) > 1 or expression.count('invT') != len(values):
                raise ValueError(f'Unsupported native activation expression at reaction {index + 1}')
            native = float(values[0]) if values else 0.0
            expected = float(energy) / rc
            if not math.isfinite(expected) or not math.isclose(native, expected, rel_tol=5.e-14, abs_tol=0.):
                raise ValueError(f'Native activation mismatch at reaction {index + 1} {branch}')
            mapping.append({'source_index': index + 1, 'native_index': index,
                            'source_equation': source_equation,
                            'native_equation': equation.strip(), 'branch': branch,
                            'source_Ea_cal_mol': float(energy),
                            'activation_temperature_K': expected,
                            'native_activation_temperature_K': native})
    if len(mapping) != 23:
        raise ValueError('Unexpected LiDryer activation branch count')
    return mapping


def prepare_source(source, header):
    if source.count('activation-energy: cal/mol') != 1:
        raise ValueError('Expected the pinned calorie source units')
    matches = re.findall(r'\bruc\s*=\s*([-+0-9.eE]+)\s*;', header)
    if len(matches) != 1 or not math.isfinite(float(matches[0])) or float(matches[0]) <= 0:
        raise ValueError('Invalid native calorie gas constant')
    rc = float(matches[0])
    activation_mapping(source, header, rc)
    temperatures = [float(v) / rc for v in re.findall(r'\bEa:\s*([-+0-9.eE]+)', source)]
    if len(temperatures) != 23:
        raise ValueError('Unexpected LiDryer activation coefficient layout')
    changed = source.replace('activation-energy: cal/mol', 'activation-energy: K')
    changed = re.sub(r'\bEa:\s*([-+0-9.eE]+)',
                     lambda m: f'Ea: {float(m[1]) / rc:.17e}', changed)
    changed = ('# Diagnostic native activation basis; generated from pinned LiDryer.\n'
               '# Canonical source remains unchanged; see activation-basis.json.\n' + changed)
    return changed, rc, temperatures


def verify_activation_only(canonical, variant, energy_ratio):
    """Refuse any change outside the source-derived activation conversion."""
    for field in ('species', 'thermo', 'transport'):
        if variant[field] != canonical[field]:
            raise ValueError(f'Unintended {field} change')
    changed = []

    def compare(before, after, path):
        if isinstance(before, dict):
            if not isinstance(after, dict) or before.keys() != after.keys():
                raise ValueError(f'Reaction layout changed at {path}')
            for key in before:
                compare(before[key], after[key], path + [key])
        elif isinstance(before, list):
            if not isinstance(after, list) or len(before) != len(after):
                raise ValueError(f'Reaction count changed at {path}')
            for index, (a, b) in enumerate(zip(before, after)):
                compare(a, b, path + [index])
        elif path[-1] == 'Ea':
            if not math.isclose(after, before * energy_ratio, rel_tol=5.e-14, abs_tol=0.):
                raise ValueError(f'Activation conversion changed at {path}')
            changed.append(path)
        elif before != after:
            raise ValueError(f'Non-activation reaction parameter changed at {path}')

    compare(canonical['reactions'], variant['reactions'], ['reactions'])
    if len(changed) != 23:
        raise ValueError('Unexpected number of converted activation energies')
    return changed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / 'references/pelec_tgreact_amr.json').read_text())
    files = {name: args.reference_source.resolve() / f'Submodules/PelePhysics/Mechanisms/LiDryer/{name}'
             for name in ('mechanism.H', 'mechanism.yaml')}
    for name, path in files.items():
        key = f'Submodules/PelePhysics/Mechanisms/LiDryer/{name}'
        if sha(path) != manifest['source_sha256'][key]:
            raise ValueError(f'Pinned source SHA mismatch: {name}')
    if sha(root / 'mechanisms/lidryer.yaml') != sha(files['mechanism.yaml']):
        raise ValueError('Canonical source differs from pinned native source')
    canonical = json.loads((root / 'mechanisms/lidryer.json').read_text())
    import cantera as ct
    if ct.__version__ != canonical['source']['runtime_cantera_version'] or ct.__git_commit__ != canonical['source']['runtime_cantera_git_commit']:
        raise ValueError('Use the canonical pinned Cantera importer environment')
    content, rc, temperatures = prepare_source(files['mechanism.yaml'].read_text(), files['mechanism.H'].read_text())
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = output / 'lidryer_native_activation.yaml'
    bundle = output / 'lidryer_native_activation.json'
    generated = output / 'lidryer_native_activation.F90'
    source.write_text(content)
    command = [sys.executable, str(root / 'tools/ingest_cantera_mechanism.py'),
               '--input', str(source), '--output', str(bundle), '--phase', 'gas']
    for key in ('module_name', 'loader_name', 'kernel_name', 'jacobian_name',
                'thermo_loader_name', 'transport_loader_name', 'symbol_prefix', 'chemistry_integrator'):
        command.extend(['--' + key.replace('_', '-'), canonical[key]])
    command.extend(['--description', 'Explicit pinned PelePhysics LiDryer activation-basis diagnostic'])
    commands = [command,
                [sys.executable, str(root / 'tools/check_mechanism_provenance.py'), '--bundle', str(bundle), '--regenerate'],
                [sys.executable, str(root / 'tools/generate_elementary_mechanism.py'), '--input', str(bundle), '--output', str(generated)]]
    record = {'qualification': False, 'status': 'preparing', 'native_R_cal_mol_K': rc,
              'activation_mapping': activation_mapping(files['mechanism.yaml'].read_text(),
                                                       files['mechanism.H'].read_text(), rc),
              'activation_temperatures_K': temperatures, 'pelec_commit': manifest['pelec_commit'],
              'inputs': {str(p): sha(p) for p in [*files.values(), root / 'mechanisms/lidryer.yaml',
                         root / 'mechanisms/lidryer.json', root / 'references/pelec_tgreact_amr.json', Path(__file__),
                         root / 'tools/ingest_cantera_mechanism.py', root / 'tools/generate_elementary_mechanism.py']},
              'commands': []}
    def save():
        (output / 'activation-basis.json').write_text(json.dumps(record, indent=2) + '\n')
    save()
    for index, cmd in enumerate(commands):
        with (output / f'command-{index}.log').open('w') as log:
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        record['commands'].append({'command': cmd, 'returncode': result.returncode})
        record['status'] = 'preparing' if result.returncode == 0 else 'failed'
        save()
        if result.returncode:
            raise SystemExit(f'Generation failed; inspect {output}/command-{index}.log')
    try:
        variant = json.loads(bundle.read_text())
        # Canonical Cantera converts cal/mol to J/kmol with exactly 4184.
        record['activation_energy_ratio'] = ct.gas_constant / (4184. * rc)
        record['verified_activation_paths'] = verify_activation_only(
            canonical, variant, record['activation_energy_ratio'])
    except (ValueError, KeyError, TypeError) as error:
        record['status'] = 'failed'
        record['error'] = str(error)
        save()
        raise
    record['status'] = 'generated'
    record['outputs'] = {str(p): sha(p) for p in (source, bundle, generated)}
    save()
    print(output)


if __name__ == '__main__':
    main()
