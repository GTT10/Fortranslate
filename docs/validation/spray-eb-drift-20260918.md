# Specular EB particle drift — 2026-09-18

Implementation: `8d96f0e3c5c42627ad2d5198b3b68f676da55cb7`, based on `cdd3011`.

`drift_spray_level_set` connects the nodal-level-set segment query to complete
parcel trajectories at constant velocity between stationary specular contacts.
It reflects normal velocity, follows the remaining flight time and permits
multiple contacts. A private population and interval inventory commit only
when every parcel succeeds. Diagnostics remain unchanged on failure. The
returned mass/momentum/energy change is the particle inventory difference;
its negative is received by the wall. Thermodynamic and composition fields
are untouched. No liquid-property table is needed for elastic reflection.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled and
selected H2/O2 tests-disabled configurations, unchanged compiler and options.
Compiler `/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_drift --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(drift|intersection)$' --output-on-failure --output-log build/spray-eb-drift-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_mechanism_provenance.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, total **0.06 s**; both individual times round to **0.00 s**.
The new drift test checks an oblique reflection with preserved tangential motion,
remaining-time position, represented impulse, elastic energy, unchanged liquid
mass/temperature, whole-population rollback after a later solid start, contact
budget rejection, two reflections in a thin fluid cavity and zero-time drift.
The existing intersection test also passed. Selected runtime compilation passed;
no selected runtime simulation was run. Provenance, Python syntax, project
contract and whitespace checks passed. No failures or skips; real-comparison
compiler warnings remain.

Artifacts: `build/spray-eb-drift-checks/{build.log,build-selected.log,ctest.log}`.

Limits: at most one million parcels and geometry cells, 1024 contacts per parcel,
and interval at most 1e6 seconds. Every proposed unreflected segment endpoint must
remain inside the root domain. Exterior boundary events are not composed with
EB events. This is a force-free serial API, not drag/heat/evaporation coupling,
liquid wall impact/film formation, distributed EB/AMR particle evolution or an
application qualification. No external parity, convergence or physical
validation was performed.
