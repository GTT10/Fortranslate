#!/usr/bin/env python3
"""Replay prepared gas families with native plane geometry and selected thermo.

Weights are reconstructed, and the EOS guess is recovered from parent held
energy; neither was recorded by QBOUND_SAMPLE. This is not a trajectory replay.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from analyze_eb_quartic_admission import read_blocks


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for flag in ('input','case','build-dir','work-dir'):
        ap.add_argument('--'+flag,type=Path,required=True)
    args=ap.parse_args();b=args.work_dir.resolve();b.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1];lib=args.build_dir.resolve()
    rows=args.case.read_text().splitlines()
    assert rows[0]=='PELEF_REACTING_EB_REFERENCE_V4'
    size=int(rows[1]);eos=int(rows[3]);dx=list(map(float,rows[5].split()));plane=float(rows[6].split()[0])
    blocks=read_blocks(args.input)
    sources=[root/'tools/reference_probes'/p for p in ('eb_extremum_slope_mod.F90','eb_bounded_tensor_mod.F90',
             'eb_quartic_kernel_mod.F90','eb_guarded_quartic_mod.F90','eb_bounded_family_replay.F90')]
    archives=[lib/'tests/liblidryer_reference_mechanism.a',lib/'libpelef_core.a']
    command=['gfortran','-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
             '-I'+str(lib/'modules'),'-I'+str(lib/'tests/lidryer_reference_bundle/modules'),
             '-I'+str(b),'-J'+str(b),*map(str,sources),*map(str,archives),'-o',str(b/'probe')]
    subprocess.run(command,check=True,capture_output=True,text=True)
    lines=[' '.join(map(str,[size,*dx,plane,eos]))]
    for block in blocks:
        assert len(block['samples'][(0,0,0)])==15
        lines.append(' '.join(map(str,block['source'])))
        for k in range(-2,3):
            for j in range(-2,3):
                for i in range(-2,3):
                    lines.append(' '.join(format(v,'.17g') for v in block['samples'][(i,j,k)]))
    stdin='\n'.join(lines)+'\n';(b/'input.dat').write_text(stdin)
    run=subprocess.run([str(b/'probe')],input=stdin,capture_output=True,text=True)
    (b/'output.dat').write_text(run.stdout);(b/'stderr.log').write_text(run.stderr)
    if run.returncode:raise RuntimeError(run.stderr)
    output=run.stdout.splitlines();assert len(output)==9*len(blocks)
    summaries=[]
    for n,block in enumerate(blocks):
        values=list(map(float,output[9*n].split()))
        assert list(map(int,values[:3]))==block['source']
        summaries.append({'source':block['source'],'theta':values[3],'mean_residual':values[4],
                          'recovered_temperature':values[5],'reconstructed_volume_fractions':values[6:]})
    report={'status':'passed','qualification':False,'families':len(blocks),'build_command':command,
            'weights_provenance':'native axis-plane geometry rebuilt from V4 case',
            'guess_provenance':'held parent EOS recovered from 900 K guess, not original stage guess',
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in [args.input,args.case,Path(__file__),*sources,*archives,b/'probe']},
            'max_mean_residual':max(r['mean_residual'] for r in summaries),
            'min_theta':min(r['theta'] for r in summaries),'results':summaries}
    (b/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('results','sha256','build_command')},indent=2))


if __name__=='__main__':main()
