# Distributed two-level MPI AMR CSV — 2026-09-17

## Change and boundary

Source parent: `f734d18d142defbd6fb4a75c3f13b340d0a834f7` on local `main`.
The implementation and test SHA-256 values are recorded below.

The MPI application writes coarse/fine CSV directly from local x slabs with
collective MPI file writes at deterministic fixed-width offsets. The buffer
holds one local x row. Empty fine ranks participate with zero write counts.
Existing files are truncated to the exact new size. Both paths and both levels'
thermodynamic states are checked collectively before opening either file.
Numerical state is read-only. File-system failures after opening may leave
partial output; this is not an atomic two-file transaction.

Final composite inventory, restriction mismatch and composition closure now
use rank-local calculations and global sum/max reductions. The invariant
thresholds are unchanged. Floating-point reduction order differs from the old
root diagnostic calculation; state evolution and CSV formatting are unchanged.
The serial writer remains independent and supplies the byte-comparison oracle.
Initial and checkpoint I/O still require a full hierarchy on root. This check
does not establish production scaling, physical validation, other distributed
output formats, or general multi-patch/multi-level AMR.

## Reproduction

Reused configured build `build/spray-mpi-migration-gnu-debug`: GNU Fortran 13.3,
OpenMPI 4.1.6, Debug with bounds/FPE checks, MPI and experimental spray enabled.
No compiler, MPI or backend configuration was changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_amr_reactive_3d pelef_amr_reactive_3d --parallel 4
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_regrid_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^mpi_amr_regrid_3d_np(4|8)$' --output-on-failure
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^regression_mpi_amr_balanced_boundary_touching_(dynamic_)?restart_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

## Evidence

Completed focused tests: **4/4 passed, 0 failures, 0 skips**. Unit np4/np8
passed in 2.78/2.39 s; static/dynamic application restart passed in
54.84/82.83 s. Builds, mechanism provenance, project contract and
`git diff --check` passed. Python AST syntax checks passed for the three
verification tools listed above; no Python source changed. No full matrix,
installed-binary audit or external solver/experimental comparison was run.

Focused unit checks cover ratios 2/3, empty fine ranks, unequal partitions,
regridding and file shrink/overwrite. Every successful coarse/fine CSV is byte
compared against the independent serial writer. Invalid fine state and
rank-inconsistent fine path reject before altering either existing file.
Application checks compare serial uninterrupted output with 4-rank stop and
2-rank resume, with cell-count balancing, wall transport and reacting elementary
corner cases. The dynamic check also moves the patch and balances again.

Logs and source hashes: `build/amr-mpi-distributed-csv-checks/`.
Application artifacts:

- `build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_restart_3d/run-pwm3knq7/`
- `build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_dynamic_restart_3d/run-_z1o_o9l/`

## Source hashes

```text
223a35acdb09bc4ce9eb0a79b6e70981c28632685cfb6fdbdf5cc15520a9442f  src/io/reactive_csv_io_3d_mod.F90
3eed6963e9542510bd7b762fff6666fddcb0d7b2a2bef21fa50a358fd73b1d60  src/parallel/mpi_amr_sparse_reactive_3d_mod.F90
6a758428f44d8ca21f3f20bb17a0f0dd3f2d6558ec7460eb7c2da04c5833ec1a  src/driver/mpi_amr_reactive_3d_application_mod.F90
a1795a4593ceb20f1834018258796afa82996108c6c1f78515eaba75dce631a1  tests/mpi/test_mpi_amr_regrid_3d.F90
```
