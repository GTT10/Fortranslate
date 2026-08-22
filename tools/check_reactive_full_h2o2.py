#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, math
from pathlib import Path

def load(path):
    with Path(path).open(newline='', encoding='utf-8') as f: return list(csv.DictReader(f))

def main():
    p=argparse.ArgumentParser(); p.add_argument('--one-d',required=True); p.add_argument('--two-d',required=True); a=p.parse_args()
    r1=load(a.one_d); r2=load(a.two_d)
    if len(r1)!=8 or len(r2)!=16: raise AssertionError((len(r1),len(r2)))
    fields=['temperature','pressure','Y_H2','Y_H','Y_O','Y_O2','Y_OH','Y_H2O','Y_HO2','Y_H2O2','Y_AR','Y_N2']
    ref={k:float(r1[0][k]) for k in fields}
    for rows in (r1,r2):
      for row in rows:
        for k in fields:
          v=float(row[k])
          if not math.isfinite(v): raise AssertionError((k,v))
          if abs(v-ref[k])>2e-10*max(1.0,abs(ref[k])): raise AssertionError((k,v,ref[k]))
        ys=[float(row[k]) for k in fields if k.startswith('Y_')]
        if min(ys)<-5e-13 or abs(sum(ys)-1)>5e-11: raise AssertionError(('composition',min(ys),sum(ys)))
    if ref['Y_HO2']<=0.0 or ref['Y_H2O2']<=0.0: raise AssertionError(('radicals',ref))
    print('full-H2/O2 1D/2D uniform-reactor parity: PASS')
    return 0
if __name__=='__main__': raise SystemExit(main())
