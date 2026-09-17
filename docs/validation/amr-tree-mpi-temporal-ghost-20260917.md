# MPI coarse-time ghost interpolation — 2026-09-17

## Source and scope

Base: `f246bd8564a617becbdbb2579d7e3b0430775248`, local `main`.
Adds partial coarse-future endpoints, agreed per-depth time fractions and optional
target-depth selection to MPI ghost filling. Sparse cell validation/fetch accepts
an optional maximum depth. The serial sampler and numerical interpolation remain
unchanged; inactive future levels require no payload allocation.

Verified source SHA-256:

- `src/parallel/mpi_amr_patch_tree_ghost_3d_mod.F90`:
  `0ace7e46362ddfc0813fff47f362f4d44a1e83a10adafb7a43382a0925d9742f`
- `src/parallel/mpi_amr_patch_tree_cells_3d_mod.F90`:
  `5cff4635f07d88a6e5c2e99283e37966ed9bd4d7cb1e7ce059c63ebf2034d431`
- `tests/unit/test_mpi_amr_patch_tree_ghost_3d.F90`:
  `23b2e6bc3de46d57a163b0e79a00bfcd13d2f6dd0a217efba2e51fdb7475d6ab`

## Completed verification

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing clean system-MPI configuration
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target test_mpi_amr_patch_tree_ghost_3d test_mpi_amr_patch_tree_cells_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_mpi_amr_patch_tree_(ghost|cells)_3d_[23]$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed and isolated selected runtime builds passed. Initial four tests passed;
the test was then extended to exercise temporal PCM/MC and zero-depth empty
endpoint queries in addition to minmod. Final CTest **4/4 passed**, 5.17 s total:
cells 2/3 ranks 0.44/0.46 s, ghosts 2/3 ranks 2.12/2.12 s. No failures/skips.
Provenance, project contract and whitespace checks passed. No Python sources
changed. Existing compiler warnings remain.

The three-level fixture compares all owned active ghost values exactly with the
serial sampler for PCM, MC and minmod, all-depth and depth 0/1/2 selection,
fractions zero/one and distinct 0.25/0.75 coarse fractions. Current and future
values differ; the finest future fields are absent. Boundary modes include
periodic, inflow/outflow, slip/no-slip/isothermal walls and prescribed wall motion.
All payload on rank zero exercises empty owners.

Different valid fractions, inconsistent endpoint presence and a missing coarse
endpoint are collectively rejected with output allocations/values unchanged.
Existing nonfinite-donor, physical-transform failure, empty-rank option rejection
and sparse-cell checks also pass. No tolerance/golden changes were made.

Logs: `build/amr-tree-mpi-temporal-checks/{build,final-build,ctest}.log`.
Binaries: `build/selected-tree-mpi-debug/tests/test_mpi_amr_patch_tree_ghost_3d`
and `test_mpi_amr_patch_tree_cells_3d` in the same directory.

## Limits

The application remains synchronous. Per-level hydro/transport stages, recursive
substeps and application subcycling controls are not connected yet. Query planning
still records repeated stencil reads; no scaling claim is made. No broad matrix,
installed-prefix audit, external solver parity or physical validation was run.
