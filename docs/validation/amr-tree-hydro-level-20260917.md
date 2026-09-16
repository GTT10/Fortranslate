# Single-depth tree hydro stepping — 2026-09-17

## Source and scope

Base: `3523a9e334e38e48d3efccd7467e870b30829295`, local `main`.
The new single-depth SSPRK2 API updates peer patches at one depth and returns
integrated stage fluxes, with coarse endpoint interpolation for both RK stages.
Inactive levels remain untouched. It is a subcycling primitive; no recursive
scheduler, interval reflux or application mode is claimed by this change.
The existing common-timestep entry point remains separate.

SHA-256 at verification:

- `src/amr/amr_patch_tree_hydro_3d_mod.F90`:
  `1a158800cd44a6118b01f055bc6e8b633d87c5dfdfab8b4cc6784c4434e946a5`
- `tests/unit/test_amr_patch_tree_hydro_level_3d.F90`:
  `910890cfe1f4c27dbfc6afd8082bf3b0cbac272441ea41e80f0ede9eeab3a4a6`

## Completed checks

GNU Fortran 13.3.0 Debug with bounds checks/backtraces/FPE traps, MPI disabled;
existing `build/selected-amr-tree-debug` build directory.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target test_amr_patch_tree_hydro_level_3d pelef_selected_amr_tree_reactive_3d_runtime -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^unit_amr_patch_tree_hydro_level_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed and selected runtime builds succeeded. Focused CTest: **1/1 passed**,
1.40 s, no failures or skips. Provenance, project contract and whitespace passed.
No Python source changed; Python syntax verification was not applicable.

The unit compares root-only PCM/periodic and PLM/physical-boundary steps against
the existing synchronous method. A two-peer refined case uses coarse endpoint
fractions 0.25/0.5, checks matching peer fluxes, demonstrates a changed solution
versus frozen donors and checks inactive state/temperature preservation. For every
active cell, the returned integrated flux divergence independently reconstructs
the accepted RK update. Backward fraction intervals and a nonphysical large step
reject without changing state, temperatures or previously returned fluxes.

Logs: `build/amr-tree-hydro-level-checks/{build,ctest}.log`.
Executable: `build/selected-amr-tree-debug/tests/test_amr_patch_tree_hydro_level_3d`.

This is focused serial API evidence, not complete subcycled evolution, MPI,
configuration-matrix qualification, external PeleC parity or physical validation.
The caller still owns timestep bounds, fine-substep accumulation, restriction,
coarse/fine flux synchronization, composite budgets and full-interval rollback.
