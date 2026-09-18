# Serial EB per-port nozzle and size models — 2026-09-18

Implementation `d9823c0`, base `a0c5132`.

The serial EB namelist adds `port_nozzles(16)` and `port_sizes(16)`. Each field
inherits its shared counterpart unless explicitly overridden. Resolved controls
feed the existing deterministic cone/size samplers and normalized restart
context. This permits coexisting solid/hollow outlets and distinct particle-size
distributions and weighting rules without changing injection, property or
checkpoint numerical kernels.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Executed tests were serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_port_models_application|spray_eb_port_profiles_application)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_eb_spray_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, total **2.18 s** (1.07 s each), no failures or skips.
The new driver performs eight launches: full, stopped, resumed and invalid-model
for each pairing. Shared solid/Weibull/mass weighting is paired with an overridden
hollow/lognormal/number-weighted second port; the reverse pairing also runs.
Unspecified swirl and diameter bounds inherit shared settings. Different port
profiles, radii, discharge coefficients and synthetic compositions remain active.

Fresh final-time records verify rim/interior geometry, swirl, per-port size
weights, component masses, integrated doses and rate-derived velocities. Full
and resumed gas/particle/film CSVs and printed budgets are byte-identical.
An invalid second-port nozzle profile rejects injection before CSV output.
The existing profile driver verifies unchanged shared-model inheritance.

Fixed/selected executable builds, provenance, Python syntax, project contract
and whitespace checks passed. Selected execution, other distribution families,
broad matrix, installed audit, external solver parity and physical nozzle
validation were not run.

Logs: `build/spray-eb-port-models-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
New fixture artifacts:
`build/spray-amr-reactive-debug/tests/spray_eb_port_models_application/nozzle__qrm165w`.
