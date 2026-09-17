# Direct MPI patch-owner migration — 2026-09-17

## Source

Base `219485daea228b0b18c8e704d5514ef447e5c243`, local `main`.
Adds same-geometry U/T migration with collective preflight and candidate publication.
No numerical kernels, mechanisms, tolerances or golden hashes changed.

SHA-256:

- `src/parallel/mpi_amr_patch_tree_distribution_3d_mod.F90`:
  `e2fb8dee7739849fc3207f1e36c6f8a69fa1d42a2fe519961a644db4d2a5592e`
- `tests/unit/test_mpi_amr_patch_tree_distribution_3d.F90`:
  `a9d06bae2eec43d3da8da90b9dddeb34b4f0e1928e0bff8cddc8fbda8559324b`

## Commands and results

GNU Fortran 13.3.0 Debug, OpenMPI 4.1.6; existing matching MPI build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_mpi_amr_patch_tree_distribution_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_mpi_amr_patch_tree_distribution_3d_[23]$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Build passed, including fixed MPI support and selected runtime. Existing compiler
warnings remain. CTest **2/2 passed**, no failures/skips, 0.94 s total (0.45/0.46 s).
Provenance, project contract and whitespace checks passed. No Python source changed.
Logs: `build/amr-tree-mpi-redistribution-checks/{build,ctest}.log`.

Three levels with refinement ratios 2 and 3 exercise split/concentrated ownership,
empty ranks, shifted owners, exact per-cell/component U and T, reverse migration,
identity migration, inconsistent new-owner maps and invalid temperature rejection.
Prior scatter/gather and partition assertions remain in the same focused tests.

## Limits

Library migration only; no application trigger, measured work policy, performance
benchmark, MPI fault recovery, full matrix, installed audit, external solver parity
or physical validation. Preflight failures preserve both caller fields; MPI transport
errors remain subject to the communicator error handler.
