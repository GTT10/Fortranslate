# Single-grid EB SDC checkpoint/restart — 2026-09-27

## Scope and source

Implementation: `31352662`, based on `437602a9` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Version remains 0.246.0.

Schemas 10/11/12 append SDC iteration count, chemistry-source presence,
external force, cumulative forcing inventory and optional species/energy
reaction estimates after gas records. The reader validates the complete file
before publishing any state, diagnostics or source arrays. Legacy split
schemas 1–6 remain separate. The single-grid fixed/selected application routes
these optional fields; the SDC checkpoint input rejection has been removed.

## Completed focused checks

- Provenance preflight passed with `mechanisms/h2o2_full.json`.
- Fixed app and checkpoint unit targets built successfully. Existing compiler
  warnings remain; this is not a warning-free or selected-binary audit.
- `unit_reactive_eb_3d_checkpoint`: 1/1 passed, 0.22 s total. Existing legacy
  checks plus schema-12 exact source round trips with chemistry on/off,
  iteration mismatch rejection and a nonfinite final reaction record. Failed
  reads preserve gas, geometry, diagnostics, exchanges and SDC sources.
- Fixed application, `cases/reactive_eb_3d/sdc.nml`: continuous 19 steps to
  2e-9 s versus stop at step 9 (9.8037685375372423e-10 s) and restart.
  Final CSV bytes and all component inventory records match exactly.
  Continuous maximum invariant error: 8.8404546598065317e-16;
  elemental error: 7.0696889804678506e-16; species closure error: zero.
- Python syntax, project contract and `git diff --check` passed.

Artifacts and the trajectory runner: `build/sdc-eb-restart-checks/`.

Commands, from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_3d test_reactive_eb_3d_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_eb_3d_checkpoint$' --output-on-failure
python3 build/sdc-eb-restart-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-eb-restart-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

No selected application run/build, EB AMR/MPI integration, full matrix,
installed-binary audit or physical validation was performed in this step.
The earlier application evidence's checkpoint rejection describes the previous
implementation and is superseded by this record.
