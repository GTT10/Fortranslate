# Synchronous MPI tree molecular transport — 2026-09-17

## Source and scope

Base: `206d3ca713e75bef29040c0b6d77e553b925f2ec`, local `main`.
Connects owner-local molecular-transport SSPRK2 to distributed viscous ghosts,
physical face kernels, reflux, state restriction and temperature recovery.

SHA-256 at verification:

- `src/parallel/mpi_amr_patch_tree_transport_3d_mod.F90`:
  `e132e566897ececbdbc7ce7428a9a158c444444f2463201c5b230f24188105d5`
- `tests/unit/test_mpi_amr_patch_tree_transport_3d.F90`:
  `f7dcdadafda392aed2cb8a0245c4b50d6ed8e9544d553ba385fcb3e1cd63cbc9`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6. Incremental build in the existing clean
system-MPI configuration; no compiler/backend option changes.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_transport_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_transport_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final build succeeded. Focused CTest **2/2 passed**, 9.16 s total (2 ranks 4.63 s,
3 ranks 4.50 s), no failures/skips. Provenance, contract and whitespace checks
passed. No repository Python source changed; syntax checks were not applicable.
An initial overlong declaration failed compilation and was wrapped before the
successful build/tests. Exact-control comparisons emit real-comparison warnings.

Fixtures use elementary H2/O2, spatially varying conserved states, three levels
with ratios 2/3, distributed/all-on-one ownership and empty ranks. All four transport
terms are active. Four modes cover PCM/MC coarse donors, periodic boundaries,
moving isothermal walls and prescribed inflow/outflow. Complete owned states,
temperatures, global boundary budget, maximum reflux correction and minimum theta
match the serial API exactly after one step. Nonowner fields remain absent.
Invalid transport diameter on one owner and rank-disagreeing dt reject collectively
and preserve prior states, temperatures and all diagnostics.

Logs: `build/amr-tree-mpi-transport-checks/{build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_transport_3d`.

## Remaining scope

The caller supplies dt and replicated species/transport records. Distributed
transport timestep bounds, full reaction/transport/hydro composition, coarse-time
endpoints, selected-MPI tree application and production scaling remain. This is
operator parity, not external solver parity or physical validation.
