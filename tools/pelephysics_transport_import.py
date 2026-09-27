"""Import cubic Simple-transport fits from an exact pinned generated header."""
import hashlib
import json
import math
from pathlib import Path
import re


def import_simple_transport(header: Path, species: list[str], manifest: dict) -> dict:
    entry = manifest['simple_transport_header']
    content = header.read_bytes()
    if hashlib.sha256(content).hexdigest() != entry['sha256']:
        raise ValueError('Simple transport header SHA mismatch')
    blob = hashlib.sha1(b'blob '+str(len(content)).encode()+b'\0'+content).hexdigest()
    if blob != entry['git_blob']:
        raise ValueError('Simple transport header Git blob mismatch')
    text = content.decode('ascii')
    for index, name in enumerate(species):
        if not re.search(r'^#define '+re.escape(name)+r'_ID\s+'+str(index)+r'\s*$', text, re.M):
            raise ValueError('Simple transport species ordering mismatch')
    n = len(species)
    if n != manifest['species_count']:
        raise ValueError('Simple transport species count mismatch')

    def coefficients(symbol, count):
        pairs = re.findall(r'\b'+symbol+r'\[(\d+)\]\s*=\s*([-+0-9.Ee]+)\s*;', text)
        if len(pairs)!=count or [int(i) for i,v in pairs]!=list(range(count)):
            raise ValueError(f'Invalid {symbol} coefficient layout')
        values = [float(v) for i,v in pairs]
        if not all(math.isfinite(v) and abs(v)<=1e6 for v in values):
            raise ValueError('Invalid Simple transport coefficient')
        return [values[i:i+4] for i in range(0,count,4)]

    return dict(model='pelephysics-simple', species=species,
                source=dict(file=entry['basename'], sha256=entry['sha256'],
                            git_blob=entry['git_blob'], commit=manifest['commit']),
                viscosity=coefficients('COFETA',4*n),
                conductivity=coefficients('COFLAM',4*n),
                binary_diffusion=coefficients('COFD',2*n*(n+1)))


def regenerate_simple_transport(bundle: dict, directory: Path) -> dict:
    """Re-read pinned bytes; never trust coefficients copied from bundle JSON."""
    manifest = json.loads((Path(__file__).resolve().parents[1]/'references/heptane_3sp.json').read_text())
    source = bundle['simple_transport']['source']
    entry = manifest['simple_transport_header']
    if source != dict(file=entry['basename'],sha256=entry['sha256'],
                      git_blob=entry['git_blob'],commit=manifest['commit']):
        raise ValueError('Unsupported Simple transport provenance')
    path = directory/source['file']
    if path.resolve().parent != directory.resolve():
        raise ValueError('Simple transport source leaves bundle directory')
    return import_simple_transport(path,bundle['species'],manifest)
