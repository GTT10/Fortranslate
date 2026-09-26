# Serial EB tree redistribution input and restart — 2026-09-27

Implementation `11e4006d`, based on `5ee9d226` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The shared serial tree application now accepts `redistribution='flux_redist'`
with common-step SDC. `none` remains the default. The selector reaches the
coupled library driver and the existing redistribution field in embedded gas
checkpoint records. Restart compares that identity before publishing restored
state. No checkpoint version change is needed.

The embedded-patch serializer now accepts flux redistribution with molecular
transport and a plane outside an individual patch. These exceptions require
explicit tree-patch context; standalone checkpoint restrictions remain. Stored
geometry validation and mechanism/configuration comparisons remain in force.

The application rejects flux redistribution with split or recursive stepping,
or with the MPI executable, because those routes are not connected yet. The
public fixture is `cases/reactive_eb_tree_3d/sdc_flux_redist.nml`. Fixed/selected
serial applications share the changed source; only the fixed build was checked.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance passed; the serial application and tree-checkpoint
  unit-test targets built successfully.
- CTest `unit_eb_tree_checkpoint_3d`: 1/1 passed, 0.52 s test / 0.56 s total.
- Three application runs completed: continuous two-step integration to 2e-10 s,
  stop/checkpoint after one step, and checkpoint continuation. The three-level,
  five-patch case enables chemistry, transport, external force and flux
  redistribution. Final state and budget CSVs are byte-identical between
  continuous and resumed runs. Maximum relative budget residual is
  7.16309624393716e-16.
- A continuation requesting `none` from the `flux_redist` checkpoint rejected
  with `Invalid tree checkpoint record` and produced no result CSV.
- Enabling subcycling with redistribution rejected at input validation with
  `Tree redistribution requires common-step SDC`.
- Python syntax, project contract and `git diff --check` passed. No build or
  unexpected runtime failure occurred; the two rejection runs were intentional.
  No skipped or interrupted runs were counted.

## Reproduction and artifacts

Artifacts: `build/sdc-eb-tree-redist-app-checks/`, including `provenance.log`,
`build.log`, `ctest.log`, `contract.log`, `run.py`, `app.log`, `summary.json`,
plus isolated inputs, logs, CSVs and checkpoint under each run directory.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target pelef_reactive_eb_tree_3d test_eb_tree_checkpoint_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_eb_tree_checkpoint_3d$' --output-on-failure
python3 build/sdc-eb-tree-redist-app-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-eb-tree-redist-app-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

This is focused application/restart evidence. No selected-mechanism build,
regrid trajectory, installed-binary audit, full matrix, increased stable CFL,
physical validation or upstream AMR redistribution parity was established.
MPI and recursive redistribution and tree state redistribution remain open.
