#!/usr/bin/env python3
"""Check isolated bounded tensor mean and signed-axis covariance."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--work-dir',type=Path,required=True)
    ap.add_argument('--continuous',action='store_true')
    args=ap.parse_args();b=args.work_dir;b.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1]
    sources=[root/'tools/reference_probes'/n for n in
             ('eb_extremum_slope_mod.F90','eb_bounded_tensor_mod.F90','eb_bounded_tensor_probe.F90')]
    if args.continuous:
        original=sources[-1]
        text=original.read_text()
        text=text.replace('bounded_child_average(q,i,j,k)', 'bounded_child_average(q,i,j,k,continuous=.true.)')
        generated=b/original.name;generated.write_text(text);sources[-1]=generated
    exe=b/'probe';cmd=['gfortran','-std=f2018','-O0','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
                     '-J'+str(b),'-I'+str(b),*map(str,sources),'-o',str(exe)]
    subprocess.run(cmd,check=True,capture_output=True,text=True)
    coords=[(x,y,z) for z in range(-2,3) for y in range(-2,3) for x in range(-2,3)]
    children=[(x,y,z) for z in range(2) for y in range(2) for x in range(2)]
    rng=random.Random(1941);cases=[];groups=[]
    for _ in range(20):
        q={idx:rng.uniform(-1,2) for idx in coords}
        group=[]
        for perm in itertools.permutations(range(3)):
            for signs in itertools.product((-1,1),repeat=3):
                def transform(idx):return tuple(signs[d]*idx[perm[d]] for d in range(3))
                group.append((len(cases),perm,signs))
                cases.append([q[transform(idx)] for idx in coords])
        groups.append(group)
    data=''.join(' '.join(format(v,'.17g') for v in q)+'\n' for q in cases)
    run=subprocess.run([str(exe)],input=data,text=True,capture_output=True,check=True)
    out=[list(map(float,l.split())) for l in run.stdout.splitlines()]
    assert len(out)==len(cases)
    maxmean=max(abs(sum(row)/8-q[62]) for q,row in zip(cases,out))
    assert maxmean<1e-13,maxmean
    maxcov=0.
    for group in groups:
        reference=next(i for i,p,s in group if p==(0,1,2) and s==(1,1,1))
        for index,perm,signs in group:
            for n,idx in enumerate(children):
                mapped=tuple(idx[perm[d]] if signs[d]>0 else 1-idx[perm[d]] for d in range(3))
                other=mapped[0]+2*mapped[1]+4*mapped[2]
                maxcov=max(maxcov,abs(out[index][n]-out[reference][other]))
    assert maxcov<1e-13,maxcov
    result={'status':'passed','qualification':False,'families':len(cases),'random_fields':20,
            'signed_axis_maps':48,'max_mean_residual':maxmean,'max_covariance_residual':maxcov,
            'build_command':cmd,'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (b/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
