#!/usr/bin/env python3
"""Audit pre-chemistry EB-tree pressure and area-weighted interface fluxes.

Input is the existing V4 probe import (two nested ratio-two octant patches).
Regular controls remove the plane only: imported states are kept, including
covered coarse restriction. They are geometry ablations, not new initializers.
Uniform controls copy one fluid cell to every record. No solver policy changes.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

FIELDS = ('mass_kg_m2_s', 'normal_momentum_Pa', 'energy_W_m2')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(prefix, n):
    with Path(str(prefix) + '_initial_pressure.csv').open() as f:
        pressure = {tuple(int(r[k]) for k in ('node', 'i', 'j', 'k')): float(r['pressure_Pa'])
                    for r in csv.DictReader(f)}
    with Path(str(prefix) + '_initial_flux.csv').open() as f:
        flux = {tuple(int(r[k]) for k in ('node', 'axis', 'i', 'j', 'k')): r
                for r in csv.DictReader(f)}
    if not pressure or not flux:
        raise ValueError("Empty audit output")
    if not all(math.isfinite(p) for p in pressure.values()):
        raise ValueError("Nonfinite pressure")
    if not all(math.isfinite(float(r[k])) for r in flux.values() for k in (*FIELDS, "aperture")):
        raise ValueError("Nonfinite flux")
    groups = {}
    for key, row in flux.items():
        node, axis, *idx = key
        a = axis - 1
        if node != 2 or idx[a] not in (0, n) or (axis == 3 and idx[a] == 0):
            continue  # z lower face is a physical wall, not a periodic interface
        if float(row['aperture']) <= 0:
            continue
        coarse = [(i + 1) // 2 for i in idx]
        coarse[a] = n if idx[a] == 0 else n // 2
        ck = (1, axis, *coarse)
        groups.setdefault(ck, []).append(row)
    records = []
    for key, fine in groups.items():
        coarse = flux[key]
        ca = float(coarse['aperture'])
        fa = sum(float(r['aperture']) for r in fine) / 4
        # Fluxes are intensive; multiply by open area before comparison.
        delta = {k: sum(float(r['aperture']) * float(r[k]) for r in fine) / 4
                    - ca * float(coarse[k]) for k in FIELDS}
        records.append({'coarse_face': key, 'open_area_fraction_difference': fa-ca, **delta})
    ghosts, valid, covered = [], [], []
    jumps, mass = [], []
    for (node, i, j, k), p in pressure.items():
        if node == 1 and all(1 <= v <= n//2 for v in (i,j,k)):
            covered.append(p - 101325.)
        elif all(1 <= v <= n for v in (i,j,k)):
            valid.append(p - 101325.)
    for key, row in flux.items():
        node, axis, *idx = key
        a = axis-1
        if node != 2 or idx[a] not in (0,n) or (axis == 3 and idx[a] == 0):
            continue
        if float(row['aperture']) <= 0:
            continue
        left = idx.copy(); right = idx.copy(); right[a] += 1
        pl = pressure[(2,*left)]; pr = pressure[(2,*right)]
        ghosts.append((pl if idx[a] == 0 else pr)-101325.)
        jumps.append(pr-pl); mass.append(float(row[FIELDS[0]]))
    def maximum(v):
        return max(map(abs,v), default=0.)
    def correlation(x,y):
        mx=sum(x)/len(x); my=sum(y)/len(y)
        xx=sum((v-mx)**2 for v in x); yy=sum((v-my)**2 for v in y)
        return sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(xx*yy) if xx*yy else None
    return {'interface_coarse_faces':len(records), 'interface_fine_faces':len(jumps),
            'max_valid_pressure_offset_Pa':maximum(valid),
            'max_covered_pressure_offset_Pa':maximum(covered),
            'max_adjacent_ghost_pressure_offset_Pa':maximum(ghosts),
            'max_interface_pressure_jump_Pa':maximum(jumps),
            'pressure_jump_mass_flux_correlation':correlation(jumps,mass),
            'max_interface_mass_flux_kg_m2_s':maximum(mass),
            'max_area_fraction_mismatch':maximum(r['open_area_fraction_difference'] for r in records),
            'max_area_weighted_flux_mismatch':{k:maximum(r[k] for r in records) for k in FIELDS},
            'faces':records}


def restrict_regular_records(lines, n):
    """Volume-average all eight fine cells for each regular covered coarse cell.

    This changes only diagnostic input; the original planar import is retained.
    Stored T/history are also averaged, but the audit uses held-energy EOS and
    never consumes reaction history.
    """
    cells = {tuple(map(int, row.split()[:4])): list(map(float, row.split()[4:]))
             for row in lines}
    result = []
    for row in lines:
        key = tuple(map(int, row.split()[:4]))
        node, i, j, k = key
        if node == 1 and max(i,j,k) <= n//2:
            children = [cells[(2,x,y,z)] for z in (2*k-1,2*k)
                        for y in (2*j-1,2*j) for x in (2*i-1,2*i)]
            value = [math.fsum(v[c] for v in children)/8 for c in range(len(children[0]))]
            row = ' '.join(map(str,key)) + ' ' + ' '.join(f'{v:.17e}' for v in value)
        result.append(row)
    return result


def summarize_stencils(prefix, require_restriction=True):
    """Enforce replay and conservation; pressure overshoot is an observation."""
    def rows(suffix):
        with Path(str(prefix)+suffix).open() as f:
            data = list(csv.DictReader(f))
        if not data or not all(math.isfinite(float(v)) for r in data for v in r.values()):
            raise ValueError('Empty or nonfinite stencil audit')
        return data
    stencils=rows('_initial_stencils.csv'); restricted=rows('_initial_restriction.csv')
    fields=('halo_scaled_residual','child_mean_scaled_residual','pressure_envelope_excess_Pa')
    result={k:max(float(r[k]) for r in stencils) for k in fields}
    for k in ('gas_scaled_residual','held_scaled_residual'):
        result[k]=max(float(r[k]) for r in restricted)
    if result['halo_scaled_residual']>1e-13 or result['child_mean_scaled_residual']>1e-13:
        raise AssertionError('Stencil replay or child-average conservation failed')
    if require_restriction and max(result['gas_scaled_residual'],result['held_scaled_residual'])>1e-13:
        raise AssertionError('Covered coarse restriction disagrees with fine fluid volumes')
    result['worst_stencil']=max(stencils,key=lambda r:float(r['pressure_envelope_excess_Pa']))
    result['restriction_checked']=require_restriction
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    lines=args.input.read_text().splitlines()
    if lines[0]!='PELEF_REACTING_EB_REFERENCE_V4':
        raise ValueError('Audit requires explicit V4 EOS policy')
    if lines[3].strip()!='1':
        raise ValueError('Audit requires converged V4 EOS')
    n=int(lines[1]); start=71
    if int(lines[70])!=2*n**3 or len(lines)!=start+2*n**3:
        raise ValueError('Invalid complete V4 import')
    # Pick a top-level fluid cell above the planar wall.
    cell=next(r.split()[4:] for r in lines[start:]
              if list(map(int,r.split()[:4]))==[1,n,n,n])
    report={'qualification':False,'definition':'time-zero held EOS; no chemistry or stepping; positive-axis SI flux; fine open area divided by four',
            'input_sha256':sha(args.input),'probe_sha256':sha(args.probe),
            'script_sha256':sha(__file__), 'cases':{},'commands':[]}
    for geometry in ('cut','regular','regular_restricted'):
        for profile in ('imported','uniform'):
            data=lines.copy()
            if geometry.startswith('regular'):data[6]='-1.0 0'
            if geometry=='regular_restricted':
                data[start:]=restrict_regular_records(data[start:],n)
            if profile=='uniform':
                data[start:]=[' '.join(r.split()[:4]+cell) for r in lines[start:]]
            inp=args.output/f'{geometry}_{profile}.dat';inp.write_text('\n'.join(data)+'\n')
            for policy in ('shared','separate'):
                name=f'{geometry}_{profile}_{policy}';prefix=args.output/name
                cmd=[str(args.probe.resolve()),str(inp.resolve()),str(prefix.resolve()),'--audit-initial']
                if policy=='separate':cmd.append('--separate-momentum')
                with (args.output/f'{name}.log').open('w') as f:
                    result=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=120)
                report['commands'].append({'argv':cmd,'returncode':result.returncode,'input_sha256':sha(inp)})
                result.check_returncode()
                summary=summarize(prefix,n)
                summary['stencil_audit']=summarize_stencils(prefix,geometry!='regular')
                if summary['max_area_fraction_mismatch']>1e-12:
                    raise AssertionError('Coarse/fine geometric areas disagree')
                if profile=='uniform':
                    if summary['max_interface_mass_flux_kg_m2_s']>1e-10:
                        raise AssertionError('Uniform stationary interface has mass flux')
                    if summary['max_area_weighted_flux_mismatch']['energy_W_m2']>1e-5:
                        raise AssertionError('Uniform stationary interface energy mismatch')
                report['cases'][name]=summary
                (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    report['status']='completed'
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:{f:v for f,v in s.items() if f!='faces'} for k,s in report['cases'].items()},indent=2))


if __name__=='__main__':
    main()
