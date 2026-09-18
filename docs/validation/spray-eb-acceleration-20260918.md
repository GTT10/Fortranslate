# Serial EB spray external acceleration — 2026-09-18

Implementation `a12b22d`, base `bf277b0`.

The allocated serial EB reactive API now accepts uniform acceleration and a
separate source inventory. This extends the regular-grid kick convention in
`spray_reactive_3d_mod` to EB fluid-volume weighting. For each half step,
delta momentum is rho*delta velocity and energy change is the dot product of
mid-kick momentum with delta velocity. Parcel changes use represented mass.
Covered cells and stationary attached films are not accelerated. Gas temperature
recovery enforces the existing 1e-9 relative internal-energy/temperature check.
Acceleration is finite and bounded by 1e6 m/s² per component; nonzero force
requires an nvar-sized ledger. Publication remains whole-step transactional.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja, reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Test executed serially; selected
library compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_impact_reactive$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, 0 failures/skips, 2.31 s total. Builds, provenance,
syntax, contract and diff checks passed. Existing splash, KH, splitting,
exterior, conservation and rollback cases remain.

The new fixture applies acceleration (1000,-2000,3000) m/s² to gas and one mixed
parcel in a cut-cell domain. Total inventory closes against boundary, wall,
chemical and force inventories at existing relative 3e-12 tolerance and
1e-30 floor. Force has nonzero momentum and zero mass/species entries. An
invalid gas solver after the first kick preserves the force sentinel, ID
watermark, particle velocity and all original gas/temperature values.

Artifacts: `build/spray-eb-force-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
The fixed-array wrapper, serial application and MPI EB paths are not connected
to this option yet. No physical or external-solver validation was run.
