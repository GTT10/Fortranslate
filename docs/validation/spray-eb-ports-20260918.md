# Serial EB nozzle ports — 2026-09-18

Implementation `4886672`, base `bcd72b3`.

The serial application now selects one of up to 16 ports for each event. Ports
have separate population, position, direction, mass, speed, cone angle, reference
diameter, liquid temperature and component fractions, with shared-input defaults.
Each port constructs its effective liquid from its own mixture state. Integrated
rate-driven velocity uses the selected port's density. Nozzle radius/model,
size distribution, profile and discharge coefficient remain shared.

Particles, their children and deposited films retain material/port identity.
Particle CSV now includes `material_id`. Restart rejects out-of-range particle
or film ports and computes injection IDs using per-event populations and positive
integrated doses. Resolved inputs enter normalized restart context; snapshots
preceding the new port fields fail strict context comparison.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Executed tests were serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_ports_application|spray_eb_rate_application|spray_eb_app_splash)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_eb_spray_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **3/3 passed**, total **3.95 s**: ports 1.06 s, shared rate 1.02 s and
legacy splash 1.82 s. No failures or skips. The new port driver performs eight
launches across solid and hollow nozzles: full, stopped, resumed and invalid
component-density cases for each. Two ports inject 8 and 4 records at 0 and
20 ns, with a zero-dose middle event. Port doses are 1e-10 and 2e-10 kg, speeds
20 and 30 m/s, and synthetic H2O/H2 fractions 0.25/0.75 and 0.75/0.25.

The final 24 records have ordered global IDs and expected port identity. Fresh
records verify distinct positions, speeds, composition, mixed-density diameters,
per-port event dose and size weights. Full and resumed gas/particle/film CSVs
and printed cumulative budgets are byte-identical. The existing shared-rate
and splash/restart drivers also passed. These fixtures do not physically
validate the illustrative mixture or demonstrate multiple-port film deposition.

Fixed/selected executable builds, provenance, Python syntax, project contract
and whitespace checks passed. Selected execution, full matrix, installed audit,
external solver parity and physical injection validation were not run.

Logs: `build/spray-eb-ports-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
Port artifacts:
`build/spray-amr-reactive-debug/tests/spray_eb_ports_application/nozzle_8rbb2ogv`.
