# Serial recursive SDC redistribution — 2026-09-27

Implementation `207c2951`, based on `09411463` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The serial recursive driver now accepts `redistribution='flux_redist'` and
accumulates each accepted selected-level transfer in an interval register.
For source level s and recipient level r, the register belongs to min(s,r).
It is applied once when that coarser level finishes its interval, after child
substeps and flux synchronization. A finer interval does not consume a transfer
belonging to a still-open coarser interval.

The driver also redistributes the integrated reflux correction. Contributions
to levels coarser than the closing interval are deferred to their registers;
other contributions and the closing register are applied immediately. Temperature
recovery and fluid-volume restriction follow. Registers and gas candidates remain
private until the whole root step succeeds. The existing `none` path is retained.

Serial application input now permits `subcycling=.true.` together with
`redistribution='flux_redist'`. Split stepping remains unsupported for this
selector. MPI recursive redistribution remains explicitly rejected pending its
interval-register connection.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`; existing compiler warnings remain.

- Mechanism provenance, serial application and full-step test target built successfully.
- `unit_eb_tree_sdc_redist_recursive_3d`: 1/1 passed,
  31.14 s test / 31.19 s total. Three levels, five patches, all transport
  operators, reactions, external force, two iterations and dt 1e-10 s.
  Level step counts are `[1,2,4]`; redistribution changes the nonuniform state.
  Relative composite closure: 5.1046165375940385e-15 with redistribution;
  5.0897607707713024e-15 for the `none` comparison.
- Rejected hydro stages preserve gas, temperature, reaction estimates, ledgers,
  limiter and caller-owned step counts.
- Three application runs completed using `sdc_flux_redist.nml` plus subcycling:
  continuous two-step integration to 2e-10 s, stop/checkpoint after one step,
  and continuation. Final state and budget CSVs are byte-identical between
  continuous and resumed runs. Maximum relative budget residual is
  9.416036854481e-16.
- A restart changing redistribution to `none` rejected without a result CSV.
  Split-mode redistribution rejected at input validation.
- Python syntax, project contract and `git diff --check` passed. No unexpected
  failures, skipped tests or interrupted runs were counted.

Artifacts: `build/sdc-eb-redist-recursive-checks/`, including `provenance.log`,
`build.log`, `ctest.log`, `detail.log`, `contract.log`, `run.py`, `app.log`,
`summary.json`, and isolated inputs/logs/CSVs/checkpoint.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target pelef_reactive_eb_tree_3d test_eb_tree_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_eb_tree_sdc_redist_recursive_3d$' --output-on-failure
python3 build/sdc-eb-redist-recursive-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-eb-redist-recursive-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

This is implementation/conservation/restart evidence for interval-closure
transfers. Temporal convergence, increased stable timestep and equivalence to
upstream AMR redistribution scheduling are not established. Intermediate reflux
admissibility checks remain. No MPI recursive run, selected build, regrid
trajectory, installed audit, full matrix or physical validation was run here.
