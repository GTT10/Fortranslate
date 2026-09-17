# Indexed tree checkpoint V2 — 2026-09-17

## Source and implementation

Base `a86f1837c6e3d154a63708081f3215522184a821`, local `main`.
Adds optional V2 checkpoint saving to serial/MPI fixed/selected applications with
`checkpoint_indexed=.true.`. Default saving remains V1; readers auto-detect both.
V2 stores checked int64 relative node offsets and a payload-end offset after
metadata/topology, using canonical fixed-width cell records. The shared header
reader verifies offsets against geometry and variable count. MPI seeks directly
to the footer and owned node ranges, avoiding the V1 root payload-index scan.
Both readers verify consumed node byte counts. Format selection is excluded from
physical context, preserving rank/policy interchange. Numerical kernels,
mechanism inputs, tolerances and golden hashes remain unchanged.

SHA-256:

- `src/io/amr_patch_tree_checkpoint_3d_mod.F90`:
  `8667fd6d267b3687ed2cd39dc749b62d55ea90431ddbba76efe36999c727bde6`
- `src/parallel/mpi_amr_patch_tree_checkpoint_3d_mod.F90`:
  `a7bac5b455fa74ece1494ddb4d1c30eb1b5a362a06cec0f49c9af51eb86f1f0e`
- `tools/check_mpi_amr_tree_restart_3d.py`:
  `08ad388b640383133714ef1f8c0a9f45ab753e942f7ff1986024f220cdfacc1f`

## Commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d test_mpi_amr_patch_tree_checkpoint_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_mpi_amr_patch_tree_checkpoint_3d_[23]|regression_mpi_amr_tree_indexed_restart_3d)$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_restart_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed/selected serial/MPI builds passed. CTest **3/3 passed**, 57.61 s total:
V1 unit 2/3 ranks 0.46/0.47 s, indexed application 56.66 s. No failures/skips.
Provenance, Python syntax, project contract and whitespace checks passed. Existing
compiler warnings remain.

The indexed application fixture covers fixed/selected reacting/transport execution,
2-to-3-rank and serial restart, one-rank continuation, exact serial/MPI V2 file
bytes and continuous/restarted CSV/budget equality. An altered first node offset
must reject before payload publication and preserve outputs. Existing context,
truncation, ledger and alias guards remain active. V1 API unit checks still run
at 2/3 ranks, including empty owners, changed assignment, serial interchange and
late temperature-corruption rollback.

Logs: `build/amr-tree-indexed-checkpoint-checks/{build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_indexed_restart_3d/`
and `mpi_tree_checkpoint_preflight_{2,3}.chk` in the test directory.

## Limits

V2 uses canonical fixed-width ASCII cell records with LF; arbitrary text reflow
belongs to the V1 path. V1-only older executables cannot read V2. The new V2
application regression is static/common-step; dynamic recursive format combinations
remain a separate qualification scope. Metadata is replicated and validation
uses temporary owner-local copies. This is not a throughput/scaling benchmark,
full matrix, installed audit, external-solver parity or physical validation.
