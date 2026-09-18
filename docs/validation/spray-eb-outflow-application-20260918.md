# Serial EB reactive/application particle outflow — 2026-09-18

Implementation `c143cdb`, base `f9146d1`.

The allocated reactive wrapper and serial namelist accept optional
`particle_outflow`. All six exterior particle faces use outflow, matching the
gas operator's existing extrapolated-outflow states in
`eb_reactive_hydro_3d_mod`. Source-half particle losses join the reactive boundary
ledger; EB exchange remains in the wall ledger. Both source halves, gas update,
IDs, computational split counts and ledgers remain transactional.

Escaped zero-mass records and components persist through the existing checkpoint
payload and ID rules, and cumulative boundary history includes particle losses.
The mode enters normalized restart context. It defaults off and explicitly
rejects simultaneous impact-film mode; other exterior particle conditions are
not exposed through this gas application.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Executed tests were serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d test_spray_eb_impact_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_outflow_application|spray_eb_impact_reactive|spray_eb_restart)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Initial fixed/selected builds failed on a 132-column Fortran line limit; the
interface and test call were wrapped. Tests launched after that failed fixed
build used stale binaries and are excluded from evidence. Final rebuilds passed.
Final CTest **3/3 passed**, total **4.54 s**: outflow application 1.11 s,
reactive fixture 2.03 s, legacy restart 1.37 s. No final failures or skips.

The reactive fixture rejects a gas solver after private first-half KH escape
and verifies incoming population, components, gas, temperature, watermark and
diagnostic sentinels. Successful KH outflow clears represented inventory and
closes composite mass/species/momentum/energy against boundary, wall and chemistry
ledgers to relative 3e-12, scale floor 1e-30.

The application injects two water parcels near the upper x face. Both escape
with zero component mass and retained IDs. Continuous and resumed gas, particle
and film CSVs and printed cumulative budgets are byte-identical. The driver
performs eight launches including expected rejections for immutable publication,
changed liquid/geometry context, truncated archive and output aliasing.

Fixed/selected executable builds, provenance, Python syntax, project contract
and whitespace checks passed. Selected runtime, impact-film/outflow combination,
broad matrix, installed audit, external parity and physical validation were not run.

Logs: `build/spray-eb-outflow-app-checks/`; final evidence in `build-final.log`,
`ctest-final.log`, `build-selected-final.log`, with initial build/test logs retained.
Application artifacts:
`build/spray-amr-reactive-debug/tests/spray_eb_outflow_application/eb_spray_7kargtdm`.
