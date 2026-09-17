# Distributed tree flux restriction and reflux — 2026-09-17

## Source and scope

Base: `95b8f781b2bb6198a1d47e13edf14cca413da550`, local `main`.
Connects area-averaged child flux transfer to parent owners and the complete
per-interval MPI reflux composition: peer synchronization, uncovered-cell
correction, physical-boundary ledger, global maximum correction and state
restriction. Optional floor depth preserves coarser endpoints.

SHA-256 at verification:

- `src/amr/amr_patch_tree_reflux_3d_mod.F90`:
  `6eae365e98d29a2a2f3a7ea69df709783b77bb1d4b981d60b64b8f29e9b8d17a`
- `src/parallel/mpi_amr_patch_tree_flux_peers_3d_mod.F90`:
  `a65038b28195802b569b657be1a94237934c855ded3b174bbd1b361dcf1ee6ac`
- `src/parallel/mpi_amr_patch_tree_flux_restriction_3d_mod.F90`:
  `9af3830faec0761308b4bd915b0ba8bc4388ba90e3db5f35fbf327c5607ff253`
- `src/parallel/mpi_amr_patch_tree_reflux_3d_mod.F90`:
  `f05b9c91b5267d0d231e96d4f5ad98fdf93c659bf15c6c98e7dd00cb356bf4da`
- `tests/unit/test_mpi_amr_patch_tree_flux_restriction_3d.F90`:
  `2e5a2dff98b61403743bc173a61d7bf150a1cdf8d01ce51178e077bfd1be5e8a`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, incremental build in the existing clean
system-MPI configuration; compiler/backend options unchanged.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_flux_restriction_3d test_mpi_amr_patch_tree_flux_peers_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_flux_(restriction|peers)_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Build succeeded. Final focused CTest **4/4 passed**, 1.87 s total, no failures/skips.
Peer tests at 2/3 ranks took 0.45/0.46 s; restriction/composed-reflux tests took
0.46/0.48 s. Provenance, contract and whitespace checks passed. No repository
Python changed; syntax checks were not applicable. Existing real-comparison and
dependent compiler warnings remain; the build is not warning-free.

Fixtures use four nodes, three levels, ratios 2/3, distributed/all-on-one owner maps,
empty ranks, periodic and physical boundaries, and spatial/level-varying two-component
manufactured fluxes. Owner flux averaging interleaved with peer synchronization
matches the serial reflux effective fluxes exactly. Parent coverage marks and
unchanged leaf marks are checked, along with idempotent matching re-restriction,
late marked-parent mismatch rollback and invalid/inconsistent child depths.

The same test executable exercises complete MPI reflux at minimum depths 0/1/2.
Corrected owned states, active effective fluxes, global boundary correction and
maximum correction match the serial API exactly. Fixtures explicitly require
nonzero state correction, and the physical case requires nonzero boundary correction.
Inactive effective slots remain absent; coarser states are compared unchanged
through the serial interval oracle. A later-level peer mismatch preserves all
prior state, effective flux outputs and scalar/vector diagnostics.

Logs: `build/amr-tree-mpi-flux-restriction-checks/{build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_flux_restriction_3d`.
The serial implementation only exports existing helpers; its numerical algorithm
and tolerance are unchanged. These are manufactured operator comparisons, not
physical flow validation or an MPI application timestep qualification.

## Remaining scope

EOS recovery, distributed RK stages, coarse-time endpoints, rebalancing and tree
application integration remain open. Reflux validates finite values, not physical
admissibility. Whole-node ownership and per-interface collective gates have no
scaling result. No selected-MPI tree application, external solver parity or physical
validation is claimed. Runtime MPI failures follow the communicator error handler.
