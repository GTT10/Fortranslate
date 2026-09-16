# Distributed two-level MPI AMR checkpoint input — 2026-09-17

Parent: `73b815834298838cbb14c49f1b6794decf1db1d0` on local `main`.

MPI restart now reads coarse/fine state directly into rank-local candidates.
A metadata-only parser retains the existing schema, species/thermo, fixed or
selected chemistry identity, transport, geometry/tagging policy, time/step and
boundary-ledger checks. It computes 64-bit byte offsets without allocating
cell fields. Each rank reads its local x rows with collective MPI file reads.
EOF/trailer, returned byte counts, line boundaries and numerical validity are
checked. Temperature consistency and coarse/fine synchronization use the same
shared checks as distributed output. Committed temperatures are recovered in
the same way as the unchanged serial reader.

Geometry, ownership, state, temperatures, time, step, initial inventories,
transport diagnostics, regrid count and boundary ledger are committed only
when all ranks pass. Restart keeps the current coarse partition and derives
fine ownership from stored geometry, supporting different rank counts and
empty fine ranks. The application no longer allocates or scatters full restart
fields on root.

Metadata remains replicated. Canonical files produced by existing serial/MPI
writers use direct reads. If fixed offsets or row parsing fail, a streaming
list-directed fallback preserves support for reformatted records, keeping only
local fields and a one-cell buffer. This fallback scans all records on each
rank and makes no scaling claim. No general multi-patch/multi-level AMR or
physical-validation claim is made.

## Reproduction

Reused `build/spray-mpi-migration-gnu-debug`, GNU Fortran 13.3 / OpenMPI 4.1.6,
Debug with bounds/FPE checks. No compiler/MPI/backend change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_amr_reactive_3d pelef_amr_reactive_3d test_mpi_amr_regrid_3d test_amr_reactive_3d_checkpoint --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(mpi_amr_regrid_3d_np(4|8)|unit_amr_reactive_3d_checkpoint|regression_mpi_amr_balanced_boundary_touching_(dynamic_)?restart_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

The 4/8-rank unit runs cover ratios 2/3, unequal partitions, empty fine ranks,
patch movement and resize. Loaded states and recovered temperatures exactly
match the independent serial reader. A malformed cell owned by the last rank
rejects collectively and preserves all checked fields, diagnostics and layout.
Files with extra trailing spaces exercise fallback and retain exact results.
The existing serial checkpoint unit remains an independent format check.
Static/dynamic physical-boundary application regressions compare serial
uninterrupted CSV with 4-rank stop / 2-rank restart, including reacting wall
transport, elementary corner cases and automatic cell-count balancing.

Logs and source hashes: `build/amr-mpi-checkpoint-input-checks/`.
Compatibility code initially failed to compile due to declaration order and
an overlong Fortran line; both were corrected. Tests that ran against old
binaries following those build failures are excluded. Final verification uses
`build-final.log` and `ctest-final.log`; earlier passes are not added to counts.

Restart artifacts (latest run directories):

- `build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_restart_3d/run-y9spe0br/`
- `build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_dynamic_restart_3d/run-cq75sft4/`

## Source SHA-256

```text
7e92dc92118767f4dc45fb819248466abb24da96ea31684a5330b89815b3953e  src/parallel/mpi_amr_checkpoint_3d_mod.F90
1b3cb4edcc21683ae2770d004e79ff3b44d5e73df2075534f161bfd8fa2483ba  src/io/amr_reactive_3d_checkpoint_mod.F90
e9421dee70586ce34c4d89bf5af35bc5000148d50945c75fc630c881d0bfdf2b  src/driver/mpi_amr_reactive_3d_application_mod.F90
be32b7ab4aaf9e5527aa8e34d29108bb46889c63e845d3da01bb7a06509d5f53  tests/mpi/test_mpi_amr_regrid_3d.F90
```

Final focused verification: **5/5 CTest tests passed, 0 failures, 0 skips**.
Build, provenance preflight, project contract, Python AST syntax checks for
three verification scripts, and diff whitespace checks passed. No full
configuration matrix, installed-binary audit, external solver parity or
physical validation was run.
