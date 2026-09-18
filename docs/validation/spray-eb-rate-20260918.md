# Serial EB injection rates and velocity — 2026-09-18

Implementation `e26590d`, base `5b5372f`.

The application integrates the existing piecewise-linear mass-rate helper for
each scheduled interval, scales the dose, and optionally derives speed from
duration, mixed injection density, outlet area and discharge coefficient. Event
mass and speed feed both particle construction and injected budgets. Zero doses
consume no IDs; restart computes expected injection IDs from integrated doses.
The namelist binds profiles, intervals and velocity controls into restart context.
Events remain batches at their start time, not continuous injection.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Executed tests were serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_rate_application|spray_eb_restart)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_eb_spray_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, total **2.44 s**: rate application 1.04 s, legacy restart
1.36 s. The rate driver runs eight launches, full/stopped/resumed/invalid-radius
for solid and hollow nozzles. Both have events at 0, 10 and 20 ns, scales 1, 0,
2, and 5 ns integration intervals. The solid profile is constant 0.002 kg/s;
the hollow profile rises from 0.001 to 0.003 kg/s over 40 ns.

Each final event injects eight records with doses 2e-11 and 2.125e-11 kg,
respectively. Fresh particle speed magnitudes match the dose-derived expression
within relative 1e-12. Diameter/weighting/component checks remain enabled. The
zero-dose middle event consumes no IDs. Full and resumed gas, particle and film
CSVs and printed cumulative budgets are byte-identical. Negative rate-driven
radius rejects initialization before CSV output.

Fixed/selected executable builds, provenance, Python syntax, contract and
whitespace checks passed. No failures or skips. Selected runtime, broad matrix,
installed audit, external solver parity and physical injection validation were
not run. Independent serial nozzle ports remain open.

Logs: `build/spray-eb-rate-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
Application artifacts:
`build/spray-amr-reactive-debug/tests/spray_eb_rate_application/nozzle_r2dnvrmo`.
