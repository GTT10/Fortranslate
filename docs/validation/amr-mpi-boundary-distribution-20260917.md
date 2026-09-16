# Sparse MPI boundary-touching distribution — 2026-09-17

Source: parent `2a72adab4e641ee2cf222ad28daa50487b333f5c` plus this change.
Logs and source hashes: `build/amr-mpi-boundary-distribution-checks/`.
Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build with bounds checks.

Sparse slab ownership now accepts valid boundary-touching patch geometry.
Scatter/gather, cell-local chemistry and average-down use the existing local
arrays and collective contracts, including ranks owning no fine cells.
Hierarchy hydro and molecular transport explicitly retain their interior-patch
restriction: fine physical ghosts and domain-surface reflux are not connected
in the distributed stepping path yet. MPI application input restrictions are
unchanged. This change does not enable boundary-touching application runs.

The MPI regression covers lower and upper corners and a full-domain patch.
Scatter/gather preserves all state and temperature values exactly. Active
local reaction matches the serial AMR result exactly after average-down.
A bad reaction interval on one rank fails collectively; unsupported hydro and
transport calls fail with neutral diagnostics. Failed calls preserve both
levels and temperatures. Existing interior-patch MPI tests remain included.
Four and eight ranks exercise zero-fine ownership and one coarse column per
rank, respectively.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^mpi_amr_reactive_3d_np(4|8)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed, no failures or skips (four ranks 72.96 s, eight
ranks 63.45 s; total 136.43 s). Mechanism provenance, project contract and
diff checks passed. No Python source changed. The build succeeded with existing
real-equality warnings. No external solver or full configuration matrix was run.

This is focused distributed data/reaction evidence, not full boundary-touching
MPI stepping, external solver parity, scaling or physical validation.
