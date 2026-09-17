# Owner-local MPI tree CSV — 2026-09-17

## Source and implementation

Base `bf6a06387e1e9115b1663220a0f47da8e5183670`, local `main`.
Final CSV uses owner-local formatting and independent MPI-IO writes at checked
serial-order byte offsets. A validation/measurement pass precedes file opening;
a second pass uses bounded chunks. Final composite inventory uses the existing
serial-order compensated small-vector traversal. Full final state/temperature
gathers are removed. Initial inventory and checkpoint gathers remain unchanged.
No numerical kernel, mechanism source, tolerance or golden hash changed.

SHA-256:

- `src/parallel/mpi_amr_patch_tree_csv_3d_mod.F90`:
  `a006007d7418848c3eb5983e54b42ff07366a633038c26f11cd8fcb0c4909b80`
- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `5f96613bd7d08547d7df20745e04dc5acbb2d31732b147bed018af2e1e6634c8`
- `tests/unit/test_mpi_amr_patch_tree_csv_3d.F90`:
  `21cbb15c2603ce0c3eed91ef0047882f924e3c43239a2c892198e8eb272e67ee`
- `tools/check_mpi_amr_tree_application_3d.py`:
  `e558a206a00c7c21f3cb2273dc47ca7348d7e1f9a7258bc24fc095df25c8a49b`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_mpi_amr_patch_tree_csv_3d pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_mpi_amr_patch_tree_csv_3d_[23]|regression_mpi_amr_tree_selected_application_3d|regression_mpi_amr_tree_weighted_adaptation_3d)$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed/selected application builds passed. CTest **4/4 passed**, 50.59 s total:
selected static application 6.90 s, weighted adaptation 42.76 s, preflight
2/3 ranks 0.44/0.45 s; no failures/skips. Provenance, Python syntax, project
contract and whitespace checks passed. Existing compiler warnings remain.

Selected static branched/root-only cases match serial CSV and budget bytes at
2/3 ranks, including empty owners. Preseeded 1.1 MB output files verify truncation
when shorter results replace older data. Fixed/selected recursive PLM adaptation
covers weighted owners, serial-byte output equality, a changed-policy 2-to-3
restart and exact checkpoint bytes. The preflight test places all payload on
the last rank and rejects an invalid final cell, different time and different
path while preserving the prior file content and length.

Logs: `build/amr-tree-mpi-csv-checks/{build,ctest}.log`.
Application artifacts remain under `build/selected-tree-mpi-debug/tests/` in
`mpi_amr_tree_selected_application_3d` and `mpi_amr_tree_weighted_adaptation_3d`.
Preflight artifacts: `mpi_tree_csv_preflight_{2,3}.csv` in that test directory.

## Limits

This is output-format/integration qualification, not a filesystem throughput or
production scaling benchmark. Formatting executes twice. Actual I/O errors after
open/truncation can leave a partial file; only preflight rejection preserves the
destination. There is no atomic rename. Initial inventory and checkpoint I/O
still gather root payloads. No full matrix, installed audit, external-solver
parity or physical validation was run.
