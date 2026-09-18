# Serial EB independent port profiles — 2026-09-18

Implementation `8048118`, base `9cbae06`.

Each active port now resolves its own profile point count, times/rates, radius
and discharge coefficient. Unspecified values inherit shared controls. Explicit
zero-point profiles retain prescribed port mass; profiled events integrate their
selected history over the supplied event interval. The common rate-velocity
switch uses selected-port density, radius and discharge coefficient. Radius also
feeds the finite outlet sampler. Resolved inputs are bound into restart context.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Executed checks were serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_port_profiles_application|spray_eb_rate_application|spray_eb_ports_application)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_eb_spray_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **3/3 passed**, total **3.25 s**: independent profiles 1.09 s, prescribed
ports 1.07 s, shared rate 1.05 s. No failures or skips. Each driver performs
eight launches across solid and hollow outlet modes, including stopped/resumed
and expected-invalid cases.

The new case inherits a constant 0.002 kg/s profile at port 1 and overrides port
2 with a linear 0.001-to-0.003 kg/s ramp over 40 ns. Radii are 0.1/0.2 mm,
discharge coefficients 0.8/0.5, and synthetic mixed-liquid compositions differ.
The final 5 ns intervals inject 1e-11 and 1.0625e-11 kg. Fresh records verify
port-dependent outlet positions, component-derived density, dose and speed;
velocity agreement is relative 1e-12. A zero-scaled middle event consumes no IDs.
Continuous and resumed CSVs and printed cumulative budgets are byte-identical.

Fixed and selected executable builds, provenance, Python syntax, project contract
and whitespace checks passed. Selected execution, broad matrix, installed audit,
external solver parity and physical profile/nozzle validation were not run.
Nozzle shape and size-distribution models remain shared across serial ports.

Logs: `build/spray-eb-port-profiles-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
New fixture artifacts:
`build/spray-amr-reactive-debug/tests/spray_eb_port_profiles_application/nozzle_06rflus5`.
