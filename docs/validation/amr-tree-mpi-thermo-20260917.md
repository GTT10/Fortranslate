# MPI tree temperature recovery and CFL — 2026-09-17

## Source and scope

Base: `02c30b0fa7d5d8499af1bd28385c3249b32fab0a`, local `main`.
Adds owner-local temperature recovery, optional single-depth recovery, and global
minimum common-step CFL selection using unchanged serial thermodynamic kernels.

SHA-256 at verification:

- `src/parallel/mpi_amr_patch_tree_thermo_3d_mod.F90`:
  `17ce054c7a004a31d0c4ef14bc456336b5d59cc0e9917bf0e7427043d71e8463`
- `tests/unit/test_mpi_amr_patch_tree_thermo_3d.F90`:
  `8ba20f929cdb91971f30f9bcc9aae90c4198075d37e7ffabae9c5577390fb29f`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6. Existing clean system-MPI configuration,
incremental build; no compiler or backend configuration change.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_thermo_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_thermo_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Build succeeded. Focused CTest **2/2 passed**, 0.99 s total (2 ranks 0.47 s,
3 ranks 0.49 s), no failures/skips. Provenance, contract and whitespace checks
passed. No repository Python changes; syntax checks were not applicable. The
control-agreement comparison emits the expected real-inequality compiler warning.

Fixtures use three levels with ratios 2/3, spatially varying elementary H2/O2
conserved states, three ownership maps, empty ranks and all-on-one ownership.
All-depth temperatures and global CFL match serial calls exactly. Selected-depth
recovery preserves other temperature fields, and conserved sources remain unchanged.
Negative density on one owner rejects recovery/CFL collectively and preserves
prior output. Inconsistent depth and CFL controls also reject without updates.

Logs: `build/amr-tree-mpi-thermo-checks/{build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_thermo_3d`.

## Limits

Thermodynamic records are caller-replicated; coefficient identity is not compared
across ranks. This is common-step CFL, not subcycle scaling. Distributed RK stages,
coarse-time endpoints, selected-MPI tree application and broader mechanism/scaling
qualification remain. No external solver parity or physical validation is claimed.
