# Coupled/application tree subcycling — 2026-09-17

## Source and behavior

Base: `a6834c484da537f7ced56ef5dd8e1275660e3b63`, local `main`.
The split API gains optional subcycling; its two transport halves use recursive
ratio-squared steps and hydro uses ratio-based steps. Chemistry retains coarse
half-interval solves on exposed cells. Fixed/selected applications expose an
optional `amr_tree_timestepping` namelist, select scaled timestep bounds and include
the mode in checkpoint identity. Regrid/checkpoint cadence and retry limits refer
to accepted coarse intervals. Defaults retain synchronous operation.

Final source SHA-256:

- Split: `12f67fc12adedbf7dfc707cda0b1cf3e22421f6c9cc5a320764d635b3ca12a10`
- Application driver: `365a804da2e0452045408becd25801eff7abcc1031ec4c954d53efa5af5b228b`
- Restart checker: `89d18e54a986dbbab1b94a65000f802bf48e8055e6aab7d45adf288cece5ae96`
- Split unit: `f2e71b8b31bd69669e8c282e8a2e6d03bb5d3131858123303b52993edcdabe1f`

Files are `src/amr/amr_patch_tree_split_3d_mod.F90`,
`src/driver/amr_patch_tree_application_3d_mod.F90`,
`tools/check_amr_tree_restart_3d.py` and `tests/unit/test_amr_patch_tree_split_3d.F90`.

## Completed checks

GNU Fortran 13.3.0 Debug, bounds checks/backtraces/FPE traps; MPI disabled.
Existing build: `build/selected-amr-tree-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d test_amr_patch_tree_split_3d -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^(unit_amr_patch_tree_split_3d|regression_amr_tree_subcycle_restart_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

- Fixed/selected application and split-test builds succeeded. A final rebuild
  after indentation-only cleanup also succeeded; no numerical edits followed
  the test build.
- Focused CTest: **2/2 passed**, 152.29 s total; split unit 51.88 s,
  application restart 100.39 s. No failures or skips.
- Provenance, project contract, Python `ast.parse` syntax and whitespace passed.
- The split unit now exercises recursive flow/transport with implicit chemistry,
  Euler/element budgets, synchronization and post-source rejection rollback.
  Its explicit synchronous case and disabled-operator equality checks remain.
- Both application frontends run the reacting, transporting moving-wall fixture
  through three scheduled geometry changes. Each stops after the first coarse
  step/regrid and resumes to byte-identical fields/budgets versus its continuous
  run. Changed integration mode, physics, schedule, corrupt/truncated snapshots,
  ledger/counter mismatch and unsafe lexical output alias are rejected while
  existing outputs remain unchanged.
- A separate byte comparison confirmed fixed/selected output equality.

Final field SHA-256:
`8efee7a6eaffbed9fa23673c3c1b277c946e9fe2094d28c464731c8dd5c671a9`.
Budget SHA-256:
`c521f1f77d3d4dac30e15fcef976b9f08d6aa31f5e30abed89c248212d28f8d7`.

Logs: `build/amr-tree-subcycle-app-checks/`.
Per-case inputs, logs, checkpoint and outputs:
`build/selected-amr-tree-debug/tests/amr_tree_subcycle_restart_3d/{fixed,selected}/`.

## Limits

This is operator-split serial subcycling. Chemistry runs over coarse half-intervals,
not after each fine hydro/transport substep. The automatic-geometry path shares
the new dispatch, but this focused application check uses prescribed regridding.
The full coarse dt is conservatively limited by the transport bound even though
transport halves span dt/2. No global convergence, full configuration matrix,
installed-binary audit, external PeleC parity, physical validation or production
scaling is claimed. Tree MPI and broader qualification remain open.
