# Cost-driven MPI tree rebalancing — 2026-09-17

## Source and scope

Base `3fad6d385aa95e8df5f38696679a121b266ce545`, local `main`.
Adds collective owner-cost planning, a minimum predicted peak-load reduction gate,
and transactional publication of migrated U/T and owners. No numerical kernels,
mechanisms, tolerances or golden hashes changed.

SHA-256:

- `src/parallel/mpi_amr_patch_tree_distribution_3d_mod.F90`:
  `d4087bf6145b2f1f03bcb8bd91a7516af13fb65ad4a13c51f90a4f41a32bea64`
- `tests/unit/test_mpi_amr_patch_tree_distribution_3d.F90`:
  `e02422c7fbbb42409e7819060f158b932393ee11627cbf4525b371e79de87114`

## Commands and completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing matching build directory.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_mpi_amr_patch_tree_distribution_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_mpi_amr_patch_tree_distribution_3d_[23]$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Build passed, including fixed MPI support and selected runtime. Compiler warnings
include intentional exact real comparisons and existing conversion warnings.
CTest **2/2 passed**, no failures/skips, 0.95 s total (0.45/0.46 s).
Provenance, project contract and whitespace checks passed. No Python source changed.
Logs: `build/amr-tree-mpi-cost-balance-checks/{build,ctest}.log`.

Synthetic patch costs [9,6,3], initially all on rank zero, yield [0,1,1] on two
ranks and [0,1,2] on three. Peak predicted load decreases from 18 to 9. A 0.75
minimum gain leaves ownership unchanged; 0.25 triggers exact U/T migration.
Tests also cover ignored invalid nonowner costs, negative owner cost rejection,
rank-dependent thresholds, migration preflight failure preserving owners/state,
no-improvement and all-zero successful no-ops. Existing transfer and partition
checks remain in the same tests.

## Limits

Costs are supplied test data, not collected runtime measurements. Application
timing/cadence, performance/scaling, full matrix, installed audit, external solver
parity and physical validation were not exercised. MPI transport failure recovery
is not provided. Load estimates omit communication/migration cost and timer noise.
