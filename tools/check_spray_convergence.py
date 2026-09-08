#!/usr/bin/env python3
"""Three-level time refinement for the closed methanol spray integration case.

This is numerical self-convergence of the selected model, not a measured spray
comparison. All material parameters, mesh, chemistry tolerances and parcel
substep controls stay fixed; only the outer maximum timestep is halved.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile

import numpy as np

from check_spray_cantera import rows, run, require


def check(exe: Path, yaml: Path, work: Path) -> None:
    histories = []
    for name, dt in [('coarse',2e-7), ('medium',1e-7), ('fine',5e-8)]:
        directory=work/name;directory.mkdir()
        args=SimpleNamespace(exe=exe,yaml=yaml,phase='gas',case='methanol-spray',backend='cvode',maximum_dt=dt)
        run(args,directory)
        histories.append(rows(directory/'run.history.csv'))
    times=np.array([r['time'] for r in histories[0]])
    values=[]
    for history in histories:
        ts=np.array([r['time'] for r in history])
        require(np.all(np.diff(ts)>0), 'Nonmonotone convergence output')
        values.append(np.array([np.interp(times,ts,[r[key] for r in history])
                                for key in ('tmin','liquid_mass')]))
    coarse_error=np.max(np.abs(values[0]-values[2]),axis=1)
    medium_error=np.max(np.abs(values[1]-values[2]),axis=1)
    scaled_coarse=float(coarse_error[0]/2200.+coarse_error[1]/1.6e-9)
    scaled_medium=float(medium_error[0]/2200.+medium_error[1]/1.6e-9)
    require(scaled_coarse>1e-9, 'Refinement test is numerically trivial')
    require(scaled_medium<scaled_coarse, 'Time refinement does not reduce the history error')
    require(scaled_coarse<.02, 'Coarse solution lies outside the integration convergence envelope')
    report={'case':'methanol-spray-time-refinement', 'dt_seconds':[2e-7,1e-7,5e-8],
            'steps':[int(h[-1]['steps']) for h in histories],
            'coarse_fine_temperature_error_K':float(coarse_error[0]),
            'medium_fine_temperature_error_K':float(medium_error[0]),
            'coarse_fine_liquid_mass_error_kg':float(coarse_error[1]),
            'medium_fine_liquid_mass_error_kg':float(medium_error[1]),
            'scaled_coarse_error':scaled_coarse,'scaled_medium_error':scaled_medium,
            'error_reduction_factor':scaled_coarse/scaled_medium,
            'note':'Model self-convergence, not physical validation or proof of second-order accuracy'}
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    print('PASS: temporal refinement reduces coupled evaporation/temperature history errors')


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,required=True)
    p.add_argument('--yaml',type=Path,required=True)
    p.add_argument('--output-dir',type=Path)
    a=p.parse_args()
    if a.output_dir:
        work=a.output_dir.resolve();work.mkdir(parents=True,exist_ok=False)
        check(a.exe.resolve(),a.yaml.resolve(),work)
    else:
        with tempfile.TemporaryDirectory(prefix='pelef-spray-refinement-') as directory:
            check(a.exe.resolve(),a.yaml.resolve(),Path(directory))


if __name__=='__main__':
    main()
