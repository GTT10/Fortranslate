#!/usr/bin/env python3
"""Check regular hydro coefficient dispatch, automatic default and invalid input."""
import argparse
from pathlib import Path
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--selected',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    source=(root/'cases/reactive_hotspot_3d/ppm_mach.nml').read_text()
    source=source.replace("riemann_solver='hllc'","riemann_solver='pelec'")
    if args.selected:
        source=source.replace("thermo_model='full_h2o2'","thermo_model='selected', composition_count=3, "
                              "composition_species='H2','O2','N2', composition_mole_fractions=.3,.15,.55")
    args.work_dir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='run-',dir=args.work_dir.resolve()))
    def run(name,value=None,scheme='characteristic_ppm',success=True):
        text=source.replace('characteristic_ppm',scheme).replace('ppm_mach.csv',name+'.csv')
        if value is not None:
            text=text.replace('&reactive_3d',f'&reactive_3d\n difmag={value},')
        path=work/(name+'.nml');path.write_text(text)
        result=subprocess.run([str(args.executable.resolve()),str(path)],cwd=work,
                              text=True,capture_output=True,timeout=90)
        log=result.stdout+result.stderr
        (work/(name+'.log')).write_text(log)
        assert (result.returncode==0)==success,(name,log)
        if success:
            coefficient=float(next(x.split(':')[1] for x in log.splitlines()
                                   if x.startswith('Artificial-viscosity coefficient:')))
            default=0.1
            assert coefficient==(default if value is None or value=='-1' else float(value))
            return (work/(name+'.csv')).read_bytes()
        assert not (work/(name+'.csv')).exists()
        assert 'Artificial-viscosity coefficient' in log
    automatic=run('automatic')
    assert run('explicit_default','.1')==automatic
    assert run('automatic_sentinel','-1')==automatic
    assert run('disabled','0')!=automatic
    run('invalid','-2',success=False)
    for scheme in ('pcm','characteristic_plm'):
        automatic=run(scheme+'_auto',scheme=scheme)
        assert run(scheme+'_default','.1',scheme)==automatic
        assert run(scheme+'_zero','0',scheme)!=automatic
        assert run(scheme+'_sentinel','-1',scheme)==automatic
        assert run(scheme+'_positive','.037',scheme)!=automatic
    print(f'PASS: coefficient input, dispatch and rejection; artifacts={work}')


if __name__=='__main__':
    main()
