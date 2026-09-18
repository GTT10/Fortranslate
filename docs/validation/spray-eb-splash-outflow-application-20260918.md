# Serial EB splash/outflow application and restart — 2026-09-18

Implementation `76bb192`, base `22baa20`.

The serial fixed/selected EB spray application now accepts `wall_impacts=.true.`
with `particle_outflow=.true.` and forwards the latter to the allocated reactive
step. Existing checkpoint fields preserve wall/exterior/chemical histories,
films, surviving sparse particle IDs and the shared watermark. No checkpoint
format or numerical model changes were needed.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Fixed application tests ran serially;
selected executable was compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(splash_outflow|outflow)_application$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, 0 failures/skips, 2.94 s total (combined case 1.79 s,
existing specular outflow case 1.11 s). Builds, provenance, syntax, contract and
diff checks passed. Existing large-local-array compiler warnings remain.

The new case uses two water parcels at x=0.001250001 m, z=0.003999999 m,
200 m/s toward the embedded wall, represented dose 1e-16 kg, and maximum step
5 ns. It generates eight splash children; four escape and four survive with
IDs 4, 5, 8, 9. Two films remain. At the 10 ns checkpoint the watermark is 10
and the surviving population is four. Resume advances to 20 ns. Gas, particle
and film CSVs and all four printed inventory lines match uninterrupted execution
byte for byte; exterior mass change is negative.

The checker also retains immutable-checkpoint, changed-input, changed-geometry,
truncated-payload and output-alias rejection. Its wet-model rejection case uses
the original interior geometry. Nine launches cover the combined case, eight
the existing specular case. Library conservation/rollback evidence is recorded
separately in `spray-eb-impact-outflow-coupling-20260918.md`.

Artifacts: `build/spray-eb-splash-outflow-app-checks/` contains `provenance.log`,
`build.log`, `ctest.log`, and `build-selected.log`. Combined fixture files are in
`build/spray-amr-reactive-debug/tests/spray_eb_splash_outflow_application/eb_spray_2sbcopmc/`.

This is serial application/restart evidence, not MPI combined-mode support,
external PeleC parity or experimental impact validation.
