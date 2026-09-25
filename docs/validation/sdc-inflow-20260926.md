# SDC specified-state and Mach inflow — 2026-09-26

Implementation `26434369`, parent `728fbe09`, local main, version 0.246.0.

SDC now accepts inflow faces, optional prescribed states/temperatures and
Mach/pressure arrays. It reuses the existing state or Mach inlet ghost filler,
enables halo CTU's inflow treatment and forwards the same inlet data to the
physical transport flux API. The shared application driver passes its resolved
inlet data. Input permits inflow alongside periodic/outflow boundaries.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `pelef_reactive_3d`, passed.
- **2 application cases passed**: a 4³ nonuniform reacting hotspot with
  molecular transport and external forcing, x- inflow and other faces outflow,
  one step to 1e-9 s with two SDC iterations. The prescribed inlet uses
  1100 K, 300 kPa and 50 m/s; the Mach inlet uses 1100 K, M=0.2 and the
  pressure-extrapolation mode. Maximum x-velocity difference is
  **0.052895538782479434 m/s**. Normalized Euler inventory residuals are
  **6.95493e-18** (state) and **3.02127e-19** (Mach), with transport theta 1.
- Provenance, Python syntax, project contract and whitespace passed. No
  failures or skips; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d --parallel 2
python3 build/sdc-inflow-checks/check.py
python3 -m py_compile tools/check_project_contract.py build/sdc-inflow-checks/check.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-inflow-checks/` contains checker, inputs, CSVs,
build/run logs, provenance and report. Selected build and CTest were not
rerun. The two inlet cases have different physical input and establish active
routing plus inventory accounting, not equivalence or physical validation.
SDC walls, AMR/MPI and full external solver comparison remain open.
