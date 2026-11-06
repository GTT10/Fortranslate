#!/usr/bin/env python3
"""Build/run the isolated volume-weighted gas-family candidate against existing EOS."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir',type=Path,required=True,help='Existing GNU fixed-core build with modules/libpelef_core.a')
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--continuous',action='store_true')
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    b=args.work_dir.resolve();b.mkdir(parents=True,exist_ok=True);lib=args.build_dir.resolve()
    sources=[root/'tools/reference_probes'/n for n in
             ['eb_extremum_slope_mod.F90','eb_bounded_tensor_mod.F90','eb_quartic_kernel_mod.F90','eb_quartic_stencil_geometry_mod.F90','eb_guarded_quartic_mod.F90',
              'eb_guarded_quartic_probe.F90']]
    if args.continuous:
        original=sources[-1]
        text=original.read_text()
        text=text.replace('bounded_correction=.true.', 'continuous_correction=.true.')
        generated=b/original.name;generated.write_text(text);sources[-1]=generated
    library=lib/'libpelef_core.a';exe=b/'guarded-probe'
    command=['gfortran','-std=f2018','-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
             '-I'+str(lib/'modules'),'-I'+str(b),'-J'+str(b),*[str(p) for p in sources],str(library),'-o',str(exe)]
    record={'status':'running','qualification':False,'production_candidate_admitted':False,'build_command':command,
            'head':subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),
            'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),*sources,library]}}
    def save():(b/'guarded-record.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    with (b/'build.log').open('w') as stream:
        result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
    record['build_returncode']=result.returncode;save()
    if result.returncode:
        record['status']='build_failed';save();raise SystemExit(result.returncode)
    record['sha256'][str(exe)]=hashlib.sha256(exe.read_bytes()).hexdigest()
    result=subprocess.run([str(exe)],capture_output=True,text=True)
    (b/'guarded.log').write_text(result.stdout+result.stderr)
    record['test_returncode']=result.returncode
    record['status']='failed' if result.returncode else 'passed_focused_candidate_checks'
    match=re.search(r'PASS guarded candidate checks=\s*(\d+)',result.stdout)
    if result.returncode==0 and not match:record['status']='missing_test_receipt'
    if match:record['checks']=int(match[1])
    record['run_command']=[str(exe)];save()
    print(result.stdout,end='')
    if record['status']!='passed_focused_candidate_checks':raise SystemExit(1)


if __name__=='__main__':main()
