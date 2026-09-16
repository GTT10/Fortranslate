# Recursive tree hydro subcycling — 2026-09-17

## Source and behavior

Base: `4d2661aca1fa2574c32db72dccfd7c19924a316b`, local `main`.
The new hydro scheduler connects per-depth RK2, ancestor endpoint interpolation,
ratio-based recursive substeps, interval reflux/restriction, temperature recovery
and the composite boundary ledger. A companion CFL helper scales local stable
steps by cumulative refinement ratios to obtain a coarse-step bound. All outputs
commit only after the full coarse interval succeeds.

SHA-256 at verification:

- `src/amr/amr_patch_tree_subcycle_3d_mod.F90`:
  `7fe179156b5e7b6157261e9dc164856bd3a01202a323aa38bd0a79ae29e8ac05`
- `tests/unit/test_amr_patch_tree_subcycle_3d.F90`:
  `31477ec7e4e8e80707f1947cf76b0aacb3836eddc307fe3de29eb394704c4bb5`

## Completed checks

GNU Fortran 13.3.0 Debug, bounds checks/backtraces/FPE traps; MPI disabled.
Existing build: `build/selected-amr-tree-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target test_amr_patch_tree_subcycle_3d pelef_selected_amr_tree_reactive_3d_runtime -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^unit_amr_patch_tree_subcycle_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed and selected runtime builds succeeded. Focused CTest: **1/1 passed**,
3.52 s, no failures or skips. Provenance, project contract and whitespace passed.
No Python source changed; Python syntax verification was not applicable.

The unit exercises root-only and branched three-level geometries, with ratios
2 then 3. PCM/periodic and PLM/physical-boundary cases check the returned 1/2/6
step counts, coarse CFL scaling, root agreement with the existing synchronous
solver, composite conservation including physical boundary flux, exact parent
restriction, positive temperatures and species closure. A rejected large step
preserves all state, temperature, budget, correction and step-count outputs.

Logs: `build/amr-tree-subcycle-checks/{build,ctest}.log`.
Executable: `build/selected-amr-tree-debug/tests/test_amr_patch_tree_subcycle_3d`.

## Remaining scope

This is a serial hydro API. The application and coupled reaction/transport path
still use their common timestep. Diffusive subcycling and its timestep bounds,
coupled/application dispatch, and tree MPI remain open. No global convergence,
performance, full configuration matrix, external PeleC parity or physical
validation is claimed. Full-tree scratch copies remain; the hydro stepper accepts
caller-chosen dt and does not enforce CFL or retry internally. Its hydro CFL helper
does not include diffusion or prescribed reservoir bounds.
