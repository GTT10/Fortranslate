# Serial EB finite nozzle and size distributions — 2026-09-18

Implementation `655d5c6`, base `0d4a0cf`.

The serial EB application forwards the existing `spray_nozzle_model` and
`spray_size_options`, plus nozzle radius, into scheduled cone injection. New
namelist controls are `nozzle_radius`, `nozzle` and `sizes`. Defaults retain point
monodisperse injection. Existing liquid composition, event dose, fluid-location
validation and transactional publication remain in force. Normalized namelist
context includes these fields; prior contexts are rejected rather than silently
restored with different nozzle semantics.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled builds. Executed checks were serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_nozzle_application|spray_eb_restart)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_nozzle_application$' --output-on-failure
python3 -m py_compile tests/application/check_eb_spray_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Legacy restart passed (1.34 s). Initial nozzle testing passed solid injection
but rejected hollow points because computed radius exceeded 1e-4 m by at most
1.1e-19 m. The checker now uses its existing relative 1e-13 radius-equality gate
for the hollow rim and strict interior bounds for solid points. No solver or
physical tolerance changed. Final nozzle CTest **1/1 passed**, test 1.03 s,
total 1.07 s. Both targeted tests therefore have passing final evidence; no skips.

The new driver runs eight launches: full, stopped, resumed and invalid-radius
cases for each solid and hollow profile. Both use eight parcels per event at
0 and 20 ns, radius 0.1 mm, half-angle 0.1 rad and swirl 0.4. Solid uses Weibull
shape 2.5 and mass weighting; hollow uses lognormal spread 2 micrometres, number
weighting and angular spread 0.02. Fresh final-time parcels verify positions,
swirl, 5–20 micrometre diameters, component masses and 1e-10 kg event dose.
Full and resumed gas/particle/film CSVs and printed budgets are byte-identical.
Negative radius rejects injection before CSV output.

Fixed and selected executable builds, provenance, Python syntax, contract and
whitespace checks passed. Selected execution, other size-distribution modes,
broad matrix, installed audit, external solver parity and physical nozzle
validation were not run. Independent serial nozzle ports remain open.

Logs: `build/spray-eb-nozzle-checks/{provenance.log,build.log,build-selected.log,ctest.log,ctest-final.log}`.
Final application artifacts:
`build/spray-amr-reactive-debug/tests/spray_eb_nozzle_application/nozzle_2shtgj2r`.
