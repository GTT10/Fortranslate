# MPI shared-face flux synchronization — 2026-09-17

## Source and scope

Base: `4811238e64320c3674db17483c2a7615c661aaa0`, local `main`.
Adds per-depth shared/periodic face reconciliation between patch owners. Only
intersecting face rectangles move. Fine-replaced marks determine priority;
equal-priority values use the serial 64-epsilon agreement rule. This is the peer
stage of reflux, not yet a complete conserved-state correction.

SHA-256 at verification:

- `src/amr/amr_patch_tree_reflux_3d_mod.F90`:
  `58bed02d6765d08c399157b522814cc227e2d8d8a17cf225c19f9ecb03d73a04`
- `src/parallel/mpi_amr_patch_tree_cells_3d_mod.F90`:
  `c10e9068a80b240143192ea9c6fa10ac870647aff660b538a9eb42250e989725`
- `src/parallel/mpi_amr_patch_tree_flux_peers_3d_mod.F90`:
  `e9bad5fac7d0df63f092c32558f055ee5d573d9fe3f68b14fed6109597aaa575`
- `tests/unit/test_mpi_amr_patch_tree_flux_peers_3d.F90`:
  `65935778fae719c7bef60ece356d0f7a546a3a2f90041e03346d8fe51dfb48b8`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6. Incremental build in the existing clean
system-MPI configuration; compiler/backend options unchanged.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_flux_peers_3d test_mpi_amr_patch_tree_cells_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_(flux_peers|cells)_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Build succeeded. Focused CTest **4/4 passed**, 1.85 s total, no failures/skips.
Cell-fetch 2/3 ranks took 0.45/0.46 s; peer-flux 2/3 ranks took 0.45/0.46 s.
Provenance, contract and whitespace checks passed. No repository Python changes;
Python syntax checking was not applicable. The build retains real-comparison,
integer-conversion/count-bound and dependent optional-array descriptor warnings.

Peer tests cover offset siblings with partial transverse overlap, remote and
same-owner pairs, empty ranks, both directions of fine-mark priority, periodic
root self-pairs on Y/Z and periodic versus physical X faces. Expected complete
arrays/masks are compared exactly, including nonoverlaps, interiors and inactive
root/child levels. Equal-priority disagreement and inconsistent active depth or
periodic geometry reject collectively with all prior fluxes/masks preserved.
Cell-fetch tests retain coverage after factoring its topology/ownership preflight
from cell-shape validation. The serial reflux code change only exports its mark
type; its numerical algorithm is unchanged.

Logs: `build/amr-tree-mpi-flux-peers-checks/{build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_flux_peers_3d`.

## Limits

Parent area averaging, state divergence correction, physical boundary budget,
EOS recovery and distributed stepping remain to be connected. Full MPI reflux,
selected-MPI tree application, external solver parity, performance/scaling and
physical validation are not established by these tests. Runtime communication
failures follow the communicator's MPI error handler.
