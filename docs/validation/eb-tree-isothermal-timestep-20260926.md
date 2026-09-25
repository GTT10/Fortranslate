# Serial tree isothermal timestep — 2026-09-26

Implementation 6fc8db0, parent 7f03bfa, local main, version 0.246.0.
Serial tree timestep selection accepts wall_temperature(6). Domain-touching
faces receive configured thermal labels; internal patch faces do not. The
existing patch isothermal heat-transfer bound is combined with hydro and
molecular bounds. Recursive stepping converts each thermal interval to the
root interval using the existing squared level-ratio transport schedule.
Missing wall temperature, invalid supplied temperature or disabled thermal
transport rejects before committing dt.

## Completed check

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing fixed tests-enabled
build/spray-amr-reactive-debug. No backend changes.
Expanded tree transport test passed **1/1**, 4.04 s (CTest total 4.09 s), no
skips in the final run. Three orientations and four/five-patch trees exercise
common and recursive estimates, positive/nonincreasing bounds and rejection
of missing temperature while preserving a sentinel dt.
For z-planar recursive cases, conduction-only baseline dt is
5.6633521555547262e-11 s and the wall-limited dt is 2.5468910442830683e-11 s.
Common-step cases in this fixture remain controlled by the existing bound.

The first test attempt enabled viscosity and species diffusion as well; those
bounds dominated every case, so the assertion requiring an actual reduction
failed. The fixture was changed to conduction-only for the timestep comparison.
No numerical tolerance or solver formula was changed to satisfy the test.
Existing heated/adiabatic/periodic transport checks remain in this test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_transport_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Build/test commands were repeated after the fixture correction. Provenance,
syntax, project contract and whitespace passed. Compiler warnings remain.
Artifacts: build/eb-tree-isothermal-timestep-checks/{build.log,build-final.log,
initial-LastTest.log,LastTest.log}.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

MPI thermal timestep selection, full reacting temperature forwarding and
application input/restart remain open. No full matrix, installed audit,
external solver parity or physical validation ran. Complete PeleC parity
remains unproven.
