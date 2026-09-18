# Fixed-population EB acceleration API — 2026-09-18

Implementation `193f379`, base `a2ddbee`.

`advance_eb_spray_reactive` forwards optional `acceleration(3)` and
`source_change(nvar)` to the allocated core. Existing callers retain their
argument order and default behavior. The wrapper keeps its fixed-population
restriction; this change adds no particle insertion/removal mode.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Test executed serially; selected
library compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_reactive$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, 0 failures/skips, 0.33 s total. Builds, provenance,
syntax, contract and diff checks passed. The fixture retains reflection,
drag/heating/evaporation, chemistry and transport. Added acceleration
(1000,-2000,3000) m/s² closes composite inventory against boundary, wall,
chemical and force ledgers at existing relative 3e-12 tolerance and 1e-30 floor.
Covered-cell state/temperature remain unchanged. Missing source output rejects
before mutation; invalid gas solver preserves source sentinel, gas and velocity.

Artifacts: `build/spray-eb-force-fixed-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
This is API and focused numerical evidence, not physical or full PeleC parity.
