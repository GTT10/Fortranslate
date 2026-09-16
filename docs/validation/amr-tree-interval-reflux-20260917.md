# Interval-scoped tree reflux — 2026-09-17

## Source and behavior

Base: `e388317951ccc9509d4c11929b6cff8fa1b2a599`, local `main`.
Reflux now accepts a minimum active depth and optionally returns effective fluxes.
Coarser fields remain unchanged, inactive flux entries may be absent, and root
boundary correction is deferred until a depth-zero interval is synchronized.
Default full-tree behavior retains the original algorithm.

SHA-256 at verification:

- `src/amr/amr_patch_tree_reflux_3d_mod.F90`:
  `0aab1bea5c57420d1ee6ff0302f57c6d875fa1214de45e2e963f72bf00de9b60`
- `tests/unit/test_amr_patch_tree_reflux_3d.F90`:
  `be03384f8d1453bb837ffe6e57629e063a09bcb2f55b7d9d28faaeebd601adef`

## Completed checks

GNU Fortran 13.3.0 Debug with bounds checks/backtraces/FPE traps; MPI disabled.
Existing build: `build/selected-amr-tree-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target test_amr_patch_tree_reflux_3d pelef_selected_amr_tree_reactive_3d_runtime -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^unit_amr_patch_tree_reflux_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed and selected runtime builds succeeded. Focused CTest: **1/1 passed**,
0.01 s, no failures or skips. Provenance, project contract and whitespace passed.
No Python source changed; Python syntax verification was not applicable.

The four-level, six-node test retains independent exposed physical-boundary and
periodic conservation checks. It now also synchronizes depths 2, 1 and 0 in turn,
feeding returned effective descendant fluxes into each longer scope. Composed
state and boundary correction agree with one full-tree correction within 2e-13.
Coarser state arrays remain exact at each partial scope, inactive flux arrays are
unallocated, and partial scopes return zero root-budget correction. Reusing fully
synchronized fluxes produces no additional correction beyond 2e-13. Invalid depth,
conflicting peer flux and nonfinite flux reject while preserving state,
diagnostics and the previous optional flux output.

Logs: `build/amr-tree-interval-reflux-checks/{build,ctest}.log`.
Executable: `build/selected-amr-tree-debug/tests/test_amr_patch_tree_reflux_3d`.

This tests flux-interval composition, not actual time-subcycled hydro. All active
fluxes in one call must cover a common interval. The recursive scheduler still
needs to select substeps, accumulate interval fluxes, recover temperature and
handle full-step rollback. This is not MPI, full-matrix, external solver parity
or physical-validation evidence.
