# Serial EB tree SDC application/restart — 2026-09-27

## Source and scope

Implementation `f44b06f9`, based on `fd726430` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The shared serial fixed/selected tree application source accepts
`time_integrator`, `sdc_iterations` and `external_forcing`. Static common-step
SDC routes the coupled tree driver, retains reaction estimates between steps,
accumulates forcing separately and restores both through V3 tree records.
The budget CSV includes the independent forcing column. MPI, subcycling and
regridding remain explicitly gated until their SDC paths are connected.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Provenance preflight passed for `mechanisms/h2o2_full.json`.
- Fixed application built successfully after correcting a namelist line-length
  error. One invocation of the previous binary rejected the new input before
  that correction; it is retained as a failed invocation, not validation.
- `cases/reactive_eb_tree_3d/sdc.nml`: three levels, five patches, chemistry and
  transport enabled, force [20000, -30000, 40000], dt 1e-10 s.
  Continuous two-step run to 2e-10 s versus stop at step 1 / 1e-10 s and
  restart: final state CSV and budget CSV are byte-identical.
- Continuous maximum composite inventory error: 6.0366240671738174e-16.
- CTest `unit_eb_tree_checkpoint_3d`: 1/1 passed, 0.52 s test / 0.57 s total.
- Python syntax, project contract and `git diff --check` passed.

Artifacts and trajectory runner: `build/sdc-eb-tree-app-checks/`.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d --parallel 2
python3 build/sdc-eb-tree-app-checks/run.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_checkpoint_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py build/sdc-eb-tree-app-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

No selected build/run, MPI build/run, regridding, recursive SDC stepping,
full configuration matrix, installed-binary audit or physical validation was
performed. The tree SDC route retains unredistributed EB rates; this small
restart case is not a cut-cell stability or temporal-convergence qualification.
