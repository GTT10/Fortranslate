# Serial EB acceleration application and restart — 2026-09-18

Implementation `bf1cc7f`, base `cf98645`.

The serial EB application exposes finite, bounded `acceleration(3)` and forwards
it through both impact and specular reactive paths. Cumulative force inventory
is separate from wall, boundary and chemical exchange and is printed when force
is enabled. The generic spray checkpoint API adds optional `external_totals`
with an `_EXT_V1` marker and a full conserved-vector row after chemistry totals.
Writing validates size/finite values; reading stages the row and publishes it
only after the complete checkpoint validates. Default zero acceleration omits
the optional record. The namelist setting participates in restart context.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Selected executable compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(force_application|restart_application)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_restart$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two distinct tests passed: force application **1/1**, 1.15 s total, and legacy
restart **1/1**, 1.41 s total. The first regex matched only the force test;
`restart_application` was not a registered name. There were no test failures or
skips. Builds, provenance, syntax, contract and diff checks passed.

The force fixture uses acceleration (0,0,1e6) m/s² and checks positive z impulse.
Full versus resumed gas, particle and film CSVs and all five printed inventory
lines (including force) are byte-identical. Existing immutable publication,
changed input/geometry, truncation and alias rejection checks remain. Each case
uses eight launches. The checkpoint carries `_EXT_V1`; legacy zero-force restart
also passes. No physical validation or external-solver comparison was run.

Artifacts: `build/spray-eb-force-app-checks/{provenance.log,build.log,ctest.log,ctest-legacy.log,build-selected.log}`;
fixture `build/spray-amr-reactive-debug/tests/spray_eb_force_application/eb_spray_5k0s2j4a/`.
MPI EB acceleration and fixed-array wrapper exposure remain pending.
