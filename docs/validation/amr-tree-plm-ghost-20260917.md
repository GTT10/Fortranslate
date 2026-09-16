# Limited tree coarse-donor interpolation — 2026-09-17

## Source and behavior

Base: `9df51707f9687c8821743a56a1dd0ac0933f13a5`, local `main`.
The tree ghost API gains optional MC/minmod reconstruction of coarse donors and
species-closure slopes. Hydro enables it with characteristic-PLM reconstruction;
PCM hydro and transport retain the existing PCM donor path. Same-resolution
peers remain exact copies. Coarse-time interpolation also supplies donor stencils.

Verified SHA-256:

- `src/amr/amr_patch_tree_ghost_3d_mod.F90`:
  `d0cbdc0fac118d826e63cc942d8e05a7f3ca5c6562f2c2aaebb2b1b92e357b7a`
- `src/amr/amr_patch_tree_hydro_3d_mod.F90`:
  `33362a7294b374f5242c53c8f5d284068d406cae497d7d4eab832c788266b710`
- `tests/unit/test_amr_patch_tree_plm_ghost_3d.F90`:
  `c7c08f088eb2ab57ea1ceac350b742e5ec53da9d5d29bfcf1e68992d2bd0b3fd`

## Completed checks

GNU Fortran 13.3.0 Debug, bounds checks/backtraces/FPE traps, MPI disabled.
Existing build: `build/selected-amr-tree-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target test_amr_patch_tree_plm_ghost_3d test_amr_patch_tree_ghost_3d test_amr_patch_tree_hydro_3d pelef_selected_amr_tree_reactive_3d_runtime -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^unit_amr_patch_tree_(plm_ghost|ghost|hydro)_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed and selected runtime builds succeeded. Focused CTest: **3/3 passed**,
8.18 s total, no failures or skips. Provenance, project contract and whitespace
checks passed. No Python source changed; Python syntax verification was not
applicable.

The new test verifies a fully three-dimensional affine field throughout a fine
patch halo, temporal reconstruction against an independent analytic value,
unchanged interior values, a complete donor subdivision's average, species
closure, invalid limiter/species-count rejection and destination rollback.
The existing ghost test retains PCM boundary/peer behavior. The hydro unit
exercises PCM/PLM with periodic and physical boundaries, root serial comparison,
nested composite conservation, species closure and failure rollback.

Logs: `build/amr-tree-plm-checks/{build,ctest}.log`.
Executable: `build/selected-amr-tree-debug/tests/test_amr_patch_tree_plm_ghost_3d`.

## Limits

Slopes use face neighbors and full 3x3x3 stencil extrema. If a stencil neighbor
itself lies in a coarser region, that sample remains PCM. A common scale retains
component bounds and species slope closure; it does not prove nonlinear EOS
admissibility. Existing reactive checks still reject inadmissible steps.
Transport donor interpolation and subcycle scheduling remain open. This is not a
global second-order convergence study, full matrix, installed-binary audit,
external PeleC parity or physical validation.
