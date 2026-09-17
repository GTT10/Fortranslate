# Single-contact EB mixed-liquid impact — 2026-09-18

Implementation `dbb6dda`, base `5180e9903570fec0db9918b7d61721fd6c566e27`.
The new adapter connects the existing nodal first-contact query to the existing
mixture impact kernel using the actual facet normal. It reuses the planar
dispatcher's correlations and publishes retained, spawned and deposited state
transactionally. The wall ledger is the inventory received by the wall.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled
and selected H2O2 tests-disabled configurations were reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_impact$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, total **0.04 s**, no skips. Selected runtime library
built successfully; no selected execution was run. Provenance, Python syntax,
project contract and whitespace checks passed. Existing dependency compiler
warnings remain in the selected build log.

The first test failed because its options left surface tension and viscosity at
their zero defaults. Explicit synthetic test values (0.072 N/m, 0.001 Pa s)
fixed the input; production checks and tolerances were unchanged.

An oblique x+y=1 plane exercises deposition, rebound and four-child dry splash
for a synthetic water/peroxide mixture. First-contact fraction is checked within
1e-14. Combined mass/momentum/energy including the wall closes within relative
1e-12 (component-wise scale floor 1e-20); liquid-component inventory closes within
absolute 1e-23 kg. The test also checks outgoing splash direction, film normal,
fluid misses, invalid wall temperature and late invalid child-ID rollback.

Artifacts: `build/spray-eb-impact-checks/{build.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
These are single-contact library checks. No EB application integration, remaining
flight, film evolution, MPI ID allocation, full matrix, external parity or
physical validation is claimed.
