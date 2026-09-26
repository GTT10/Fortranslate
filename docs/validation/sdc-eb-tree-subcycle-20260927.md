# Recursive serial EB tree SDC — 2026-09-27

## Source and behavior

Implementation `e21ae53a`, based on `dc83a4fa` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`advance_reactive_eb_tree_sdc_subcycle_3d` advances parent levels first and
stores their gas and source predictor endpoints. Finer levels take refinement-
ratio substeps, interpolate ancestor data and accumulate their integrated
fluxes. Adjacent levels synchronize at matching interval endpoints. Boundary
exchange comes from the final synchronized root flux; wall, chemistry and
external ledgers sum the level-interval contributions. Tree, reaction fields,
ledgers, limiter and optional int64 level counts publish only on success.

The serial application now routes SDC with `subcycling=true`. Automatic
interval estimation has an optional coupled-subcycling policy: diffusion
bounds scale by the refinement ratio rather than its square, matching SDC's
coupled hydro/transport substep count. Legacy split policy is unchanged.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- Full-tree, transport-test and application targets built successfully.
- `unit_eb_tree_full_3d`: 1/1 passed, 75.39 s test / 75.44 s total.
  Three levels/five patches, two SDC iterations, chemistry, all molecular
  transport flags and force [20000, -30000, 40000], dt 1e-10 s.
  Counts are exactly [1, 2, 4]. Composite ledger residual is
  5.0897607707713024e-15, below the existing 2e-11 bound.
  Invalid hydro viscosity preserves state, reaction estimates, counts and
  all diagnostics after private old-transport evaluation.
- Application: continuous two coarse steps to 2e-10 s versus checkpoint at
  step 1 / 1e-10 s and restart. State and budget CSVs are byte-identical.
  Continuous composite inventory error: 6.0366203415084809e-16.
- After adding coupled automatic-interval scaling, rebuilt the application
  and transport test. `unit_eb_tree_transport_3d`: 1/1 passed, 5.21 s test /
  5.26 s total. New cases check the coupled bound lies between common and
  legacy split recursive bounds. The fixed-interval application trajectory
  above does not exercise automatic interval selection.
- Python syntax, project contract and `git diff --check` passed.

Artifacts and trajectory runner: `build/sdc-eb-tree-subcycle-checks/`.
`test-details.log` retains the recursive conservation residual.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_full_3d pelef_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_full_3d$' --output-on-failure
python3 build/sdc-eb-tree-subcycle-checks/run.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_transport_3d pelef_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_transport_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py build/sdc-eb-tree-subcycle-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

No recursive regrid trajectory, selected build, MPI run, full matrix,
installed-binary audit or physical validation was performed. Ratios other
than two, temporal convergence and large-tree scaling were not newly tested.
The tree continues to use unredistributed EB rates.
