# Sparse MPI tree cell requests — 2026-09-17

## Source and scope

Base: `93ed1eed853c57fac9f47a1ed54a3a60f2bd68e1`, local `main`.
Adds owner-directed cell requests without assembling global fields. Requests
retain caller order and duplicates, and may be empty on any/all ranks.

SHA-256 at verification:

- `src/parallel/mpi_amr_patch_tree_cells_3d_mod.F90`:
  `42f21a1f0ae240066be2d8e94d50a4b803a2c3de2464484d1899e85a22c1e3e1`
- `tests/unit/test_mpi_amr_patch_tree_cells_3d.F90`:
  `c3dc10a4daa829894291978f892050d6c7ec70ab2d0dd5f3cd445db8119a6b3d`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, incremental build in the clean system-MPI
configuration documented by [distribution evidence](amr-tree-mpi-distribution-20260917.md).

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_cells_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_cells_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Build succeeded; the integer division count-bound expression emits a constant
truncation warning. Focused CTest **2/2 passed**, 0.92 s total (2 ranks 0.44 s,
3 ranks 0.46 s), no failures or skips. Provenance, contract and whitespace checks
passed. No repository Python source changed; syntax checking was not applicable.

Coverage uses three levels, ratios 2/3 and two spatially varying components:
local/remote unordered requests, duplicates, zero-request ranks, all-empty queries,
distributed/all-on-one ownership, exact output ordering and unchanged source fields.
Invalid node/cell requests, inconsistent valid topology/owners/component counts and
nonfinite requested values reject collectively and preserve prior output.

Logs: `build/amr-tree-mpi-cells-checks/{build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_cells_3d`.

## Limits

This is a cell exchange primitive, not completed MPI ghost filling or evolution.
Donor selection, boundary transforms, coarse space/time interpolation and
restriction/reflux still require integration. Only requested data are checked for
finiteness. MPI runtime failures follow the communicator error handler. No scaling,
selected-MPI application, external parity or physical validation claim is made.
