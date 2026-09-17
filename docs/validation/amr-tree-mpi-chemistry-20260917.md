# MPI tree chemistry and reaction/hydro splitting — 2026-09-17

## Source and scope

Base: `866d4582433204083ca3f8306bc55cb4b649f451`, local `main`.
Adds exposed-cell chemistry on patch owners, followed by distributed restriction
and temperature recovery, plus transactional reaction/hydro Strang composition.

SHA-256 at verification:

- `src/parallel/mpi_amr_patch_tree_chemistry_3d_mod.F90`:
  `a32eb32f7b7637bbd8a3c432b2bedade57dd297366055c45254a27d9b9eebff5`
- `tests/unit/test_mpi_amr_patch_tree_chemistry_3d.F90`:
  `b00632dce35e491d1450f5c142842d2727c1ec966f869415f7293bd7489b04fd`

## Verification commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6; incremental build in the existing clean
system-MPI configuration, with no compiler or optional-backend changes.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_chemistry_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_chemistry_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Build succeeded. Focused CTest **2/2 passed**, 152.29 s total (2 ranks 76.21 s,
3 ranks 76.04 s), no failures/skips. Provenance, project contract and whitespace
checks passed. No repository Python source changed; syntax checks were not
applicable. Exact interval/tolerance agreement emits a real-comparison warning.

Fixtures use elementary H2/O2 with nonzero radical seeding, three levels with
ratios 2/3, distributed/all-on-one ownership, empty ranks, explicit/implicit ODE
integration, PCM/PLM hydro, periodic and physical inflow/outflow/slip-wall faces.
Standalone chemistry and R/2-H-R/2 state/temperature arrays are compared exactly
to serial APIs; leaf species changes are required to establish active chemistry.
The split budget and maximum reflux diagnostic also match serial. A hydro failure
after successful chemistry, invalid reaction on one owner, and rank-disagreeing
interval check collective rejection and state/temperature/diagnostic preservation.

Logs: `build/amr-tree-mpi-chemistry-checks/{build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_chemistry_3d`.

## Remaining scope

Full thermodynamic/reaction coefficient identity is a caller contract; the API
compares counts and integration controls, not every record. Molecular transport,
coarse-time endpoints/subcycling, selected-MPI tree application and broader
mechanism/scaling qualification remain. No external solver parity or physical
validation is established by these operator comparisons.
