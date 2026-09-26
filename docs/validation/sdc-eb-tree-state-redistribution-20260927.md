# Serial common-step tree StateRedist SDC — 2026-09-27

Implementation `5174ce00`, based on `64255828` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`redistribute_sdc_eb_tree_state_rate_3d` converts a composite state redistribution
into an SDC source: `[S(U + dt*R) - U]/dt`. Covered cells are excluded before
arithmetic, and the candidate rate is published only after successful validation.
The serial common-step SDC driver now applies this convention to hydro and
old/new molecular-transport rates, following the existing single-grid driver.
It applies the geometry-dependent linear state operator directly to the final
reflux correction, then recovers temperature and restricts the result.

Serial application input now accepts `redistribution='state_redist'` and
`state_redist_target_volume_fraction` (default 0.5). Both reach the existing
checkpoint configuration identity. Embedded tree records allow a plane outside
an individual patch for this mode; standalone checkpoint rules remain.
`cases/reactive_eb_tree_3d/sdc_state_redist.nml` provides the reacting fixture.
MPI and recursive state-redistribution application modes remain explicitly
rejected. Common/recursive FluxRedist paths remain separate.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`; existing compiler warnings remain.

- Mechanism provenance and the serial application, full-step and checkpoint
  unit targets built successfully.
- CTest: 2/2 passed, 13.98 s total. `unit_eb_tree_checkpoint_3d` 0.52 s;
  `unit_eb_tree_sdc_state_redist_3d` 13.41 s.
- Coupled library test: three levels, five patches, reactions, all transport
  operators, force, two iterations, dt 1e-10 s. State redistribution changes
  the nonuniform state. Relative composite closure is 3.4340988047347421e-15
  (`none` comparison: 3.7650552248896492e-15). A rejected hydro stage preserves
  gas, temperature, reaction estimates, inventories and limiter.
- Three application runs completed: continuous two-step integration to
  2e-10 s, checkpoint after one step, and continuation. Final state and budget
  CSVs are byte-identical. Maximum relative budget residual is
  6.036624067173818e-16.
- Three intentional rejections passed: restart changing the redistribution
  method, restart changing target fraction from 0.5 to 0.75, and unsupported
  recursive state redistribution.
- Python syntax, project contract and `git diff --check` passed. No unexpected
  failures, skipped tests or interrupted runs were counted.

Artifacts: `build/sdc-eb-tree-state-redist-checks/`, including `provenance.log`,
`build.log`, `ctest.log`, `detail.log`, `contract.log`, `run.py`, `app.log`,
`summary.json`, and isolated inputs/logs/CSVs/checkpoint.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_eb_tree_full_3d pelef_reactive_eb_tree_3d test_eb_tree_checkpoint_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_eb_tree_(sdc_state_redist|checkpoint)_3d$' --output-on-failure
python3 build/sdc-eb-tree-state-redist-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-eb-tree-state-redist-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

This is bounded coupling/conservation/restart evidence. It does not establish
upstream AMR StateRedist parity, a larger stable CFL, temporal convergence or
physical accuracy. MPI/recursive StateRedist, selected builds, regrid trajectories,
installed audits and a full configuration matrix were not checked here.
Intermediate reflux admissibility checks remain in place.
