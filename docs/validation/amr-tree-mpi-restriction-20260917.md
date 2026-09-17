# MPI tree restriction — 2026-09-17

## Source and scope

Base: `bb3bddd6b3010717c5b5b31c884d8d9229d872dc`, local `main`.
Adds owner-to-parent bottom-up averaging with collective validation and whole-call
rollback. Optional minimum parent depth supports interval-scoped synchronization.
The existing cell-fetch metadata preflight is shared with this operation.

SHA-256 at verification:

- `src/parallel/mpi_amr_patch_tree_cells_3d_mod.F90`:
  `70d96f253f9eb7144d612bed0267944fb44675c1d35f1675423b8693c154de46`
- `src/parallel/mpi_amr_patch_tree_restriction_3d_mod.F90`:
  `0b4204f9ffdda1d06e7cf1de842cb44035d2327477b6325b82d7cb816ff7a382`
- `tests/unit/test_mpi_amr_patch_tree_restriction_3d.F90`:
  `0586a66eff5afc3663b4f6c8cf03ced4421a01a253044d1ed8110403b65f86b6`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6. Incremental build in the existing clean
system-MPI configuration; no compiler or optional-backend configuration changes.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_restriction_3d test_mpi_amr_patch_tree_cells_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_(restriction|cells)_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Build succeeded. Focused CTest **4/4 passed**, 1.84 s total; no failures/skips.
Cell-fetch 2/3-rank tests took 0.45/0.45 s, restriction 0.45/0.46 s.
Provenance, project contract and whitespace checks passed. No repository Python
source changes; syntax checks were not applicable. The shared preflight retains
the existing integer-division truncation warning for its metadata count bound.

Restriction fixtures use four nodes, three levels, ratios 2/3 and two components.
They cover three owner maps (including empty ranks, all-on-one ownership and
cross-rank ancestor chains) and minimum depths 0/1/2. Every owned result matches
serial restriction exactly; unowned slots remain absent. The serial oracle's
exposed composite inventory is unchanged. Comparisons cover uncovered cells and
coarser parents as well as the averaged boxes. Mismatching/out-of-range depth,
inconsistent owners/valid topology, missing child allocation and NaN child data
reject collectively; prior fields remain unchanged. Cell-fetch regression tests
cover the extracted preflight and retain ordering/empty-request/rollback checks.

Logs: `build/amr-tree-mpi-restriction-checks/{build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_restriction_3d`.

## Limits

Only restricted boxes move between child/parent owners. This does not implement
flux reflux or distributed evolution, EOS recovery, rebalancing or tree application
integration. Per-node collective synchronization and whole-node ownership have no
scaling qualification. MPI runtime failures remain subject to the communicator's
error handler. No external solver parity, selected-MPI tree application or physical
validation is established by these focused tests.
