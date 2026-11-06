"""Read Fuego inverse bounds and liquid reference temperature from pinned bytes."""
import hashlib
import math
from pathlib import Path
import re

from pelephysics_transport_import import regenerate_simple_transport


def reference_eos_options(bundle: dict, directory: Path) -> dict:
    """Preserve header SHA/blob/species checks and the original thermo fits."""
    if bundle.get('simple_transport') != regenerate_simple_transport(bundle, directory):
        raise ValueError('Pinned Simple transport fits required')
    source = bundle['simple_transport']['source']
    data = (directory / source['file']).read_bytes()
    blob = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
    if hashlib.sha256(data).hexdigest() != source['sha256'] or blob != source['git_blob']:
        raise ValueError('Pinned EOS header SHA/blob mismatch')
    text = data.decode('ascii')
    bounds = []
    for symbol in ('GET_T_GIVEN_EY', 'GET_T_GIVEN_HY'):
        functions = re.findall(r'\b' + symbol + r'\s*\([^{}]*\)\s*\{(.*?)^\}', text, re.S | re.M)
        if len(functions) != 1:
            raise ValueError(f'Expected one pinned {symbol} definition')
        values = []
        for name in ('tmin', 'tmax'):
            literals = re.findall(r'\b' + name + r'\s*=\s*([-+0-9.eE]+)\s*;', functions[0])
            if len(literals) != 1:
                raise ValueError(f'Expected one literal {symbol} {name}')
            values.append(float(literals[0]))
        if not all(math.isfinite(value) for value in values) or not 0 < values[0] < values[1]:
            raise ValueError(f'Invalid pinned {symbol} temperature bounds')
        bounds.append(values)
    if bounds[0] != bounds[1]:
        raise ValueError('Pinned internal-energy and enthalpy inverse bounds disagree')
    return dict(pinned_eos_extrapolation=True,
                pinned_eos_temperature_min=bounds[0][0],
                pinned_eos_temperature_max=bounds[0][1],
                source=dict(source), definitions=['GET_T_GIVEN_EY', 'GET_T_GIVEN_HY'])


def reference_fuel_temperature(source: bytes) -> float:
    """The caller has already checked the complete case-input SHA."""
    values = re.findall(r'^\s*particles\.fuel_ref_temp\s*=\s*([^\n#]+)', source.decode('ascii'), re.M)
    if len(values) != 1 or len(values[0].split()) != 1:
        raise ValueError('Expected one pinned fuel reference temperature')
    value = float(values[0])
    if not math.isfinite(value) or value <= 0:
        raise ValueError('Invalid pinned fuel reference temperature')
    return value
