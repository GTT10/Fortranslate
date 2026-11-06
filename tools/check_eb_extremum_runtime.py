#!/usr/bin/env python3
"""Replay isolated scalar EP slopes on recorded EB donor lines.

No geometry, EOS, temporal evolution, or fixed-reference qualification.
"""
import argparse
import collections
import hashlib
import json
import math
from pathlib import Path
import subprocess
from analyze_eb_quartic_admission import read_blocks
from check_eb_extremum_slope import check


def mc(q):
    left, right = q[2]-q[1], q[3]-q[2]
    centered = (q[3]-q[1])/2
    return (math.copysign(min(abs(centered),2*abs(left),2*abs(right)),centered)
            if left*right > 0 else 0.)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--work-dir', type=Path, required=True)
    ap.add_argument('--bounded-quartic', action='store_true')
    args = ap.parse_args()
    build = check(args.work_dir, args.bounded_quartic)
    blocks = read_blocks(args.input)
    cases = []
    occurrences = collections.Counter()
    for block in blocks:
        source = tuple(block['source'])
        occurrences[source] += 1
        samples = block['samples']
        for component in range(len(samples[(0,0,0)])):
            for axis in range(3):
                for j in range(-2,3):
                    for i in range(-2,3):
                        indices = [(v,i,j) if axis==0 else (i,v,j) if axis==1 else (i,j,v)
                                   for v in range(-2,3)]
                        q = [samples[index][component] for index in indices]
                        cases.append({'source':source,'fill':occurrences[source],
                                      'component':component+1,'axis':axis+1,
                                      'offset':[i,j],'q':q})
    stdin = ''.join(' '.join(format(v,'.17g') for v in c['q'])+'\n' for c in cases)
    run = subprocess.run([str(args.work_dir/'probe')],input=stdin,text=True,capture_output=True,check=True)
    rows = [list(map(float,line.split())) for line in run.stdout.splitlines()]
    assert len(rows)==len(cases)
    changed = collections.Counter()
    central = []
    for c,row in zip(cases,rows):
        q = c['q']; slope,low,high = row
        tol = 2e-14*max(1.,*map(abs,q))
        assert all(map(math.isfinite,row)) and abs((low+high)/2-q[2])<=tol
        delta = slope-mc(q)
        if delta != 0:
            changed[(c['fill'],c['component'],c['axis'])] += 1
        if c['offset']==[0,0]:
            central.append({**c,'ep_slope':slope,'mc_slope':mc(q),'children':[low,high]})
    # A declared analytic comparison for the saved 16-cell shear fixture only.
    profile = []
    h = 2*math.pi/16
    for r in central:
        if r['fill']!=1 or r['axis']!=2 or r['component']!=2:
            continue
        src = r['source']
        density = next(x for x in central if x['source']==src and x['fill']==1
                       and x['axis']==2 and x['component']==1)
        assert max(density['q'])-min(density['q'])<1e-13
        rho = density['q'][2]
        expected = [10*math.sin((src[1]+k-.5)*h)*math.sin(h/2)/(h/2)
                    for k in range(-2,3)]
        assert max(abs(a/rho-b) for a,b in zip(r['q'],expected))<1e-10
        for parity in (0,1):
            fine = 2*src[1]-1+parity
            exact = 10*math.sin((fine-.5)*h/2)*math.sin(h/4)/(h/4)
            ep = r['children'][parity]/rho
            old = (r['q'][2]+(-.25 if parity==0 else .25)*r['mc_slope'])/rho
            profile.append({'source':src,'fine_y':fine,'ep_error':ep-exact,'mc_error':old-exact})
    result = {'status':'passed','qualification':False,'scope':'scalar saved-stencil replay only',
              'variant':build['variant'],
              'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
              'binary_sha256':build['binary_sha256'],
              'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'blocks':len(blocks),'lines':len(cases),
              'changed_lines':[{'fill':f,'component':c,'axis':a,'count':n}
                               for (f,c,a),n in sorted(changed.items())],
              'initial_velocity_profile':profile,'central_lines':central}
    (args.work_dir/'runtime-report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('central_lines','initial_velocity_profile')},indent=2))
    print('initial velocity max errors:',max(abs(x['ep_error']) for x in profile),max(abs(x['mc_error']) for x in profile))


if __name__=='__main__':
    main()
