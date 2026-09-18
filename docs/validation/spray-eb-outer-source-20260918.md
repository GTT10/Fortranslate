# Serial EB source exterior inventories — 2026-09-18

Implementation `60fcdf1`, base `9a4bbd3`.

The allocated source API adds the optional pair `outer_faces(6)` and
`outer_change(nvar)`. Its existing five-entry output remains the EB wall ledger.
After source exchange and optional KH release, it computes liquid internal
energy from the trial component state and invokes combined EB/exterior drift.
Escaped records retain IDs with zero mass and zero component inventories.
Mapped vapor-species losses join exterior mass/momentum/energy accounting.

KH eligibility clips the checked parent trajectory at the first exterior
crossing, so an out-of-domain endpoint no longer prevents KH plus outflow.
Source-first breakup/flight ordering remains unchanged. Gas, temperature,
parcels, component masses, IDs and both ledgers publish only after the full
interval succeeds. Simultaneous impact-film arguments are explicitly rejected;
this increment connects specular EB source flight, not film/exterior events.
Reactive/application/MPI connections remain subsequent work.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Executed test was serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact_coupling --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_impact_coupling$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, test 0.18 s, total 0.22 s. Existing deposition, splash,
film evolution, KH release and rollback assertions remain. Added cases send a
mixed parcel toward the upper x face at 1000 m/s over 1e-8 s, with source drag,
heating and evaporation disabled. With KH off, the parent escapes; with KH on,
the parent and four children escape and retain watermark 104 above incoming 100.
Both cases clear all represented component inventory and record 1e-10 kg total
loss, split 0.25/0.75 across mapped species. Composite gas/particle change equals
EB plus exterior inventories to relative 3e-12, scale floor 1e-30.

A longer interval with a one-substep limit fails after private escape and
preserves incoming gas, temperature, particle/component state, watermark and
both ledger sentinels. No failures or skips. Fixed/selected library builds,
provenance, Python syntax, contract and whitespace checks passed. Selected
runtime, broad matrix, installed audit, external parity and physical validation
were not run.

Artifacts: `build/spray-eb-outer-source-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
