# Distributed two-level MPI AMR checkpoint output — 2026-09-17

Parent: `21d18b6e02ff8268e367b08922541b7b12117aa7`, local `main`.

Checkpoint output now writes state/temperature rows directly from local slabs
with collective MPI file I/O. Metadata is serialized through the same helper
as the serial writer, preserving existing schemas, fixed/selected chemistry
identities, transport controls, dynamic geometry and physical boundary budgets.
Metadata generation uses scratch formatted files; it contains no cell payload.
The output buffer holds one local x row. Empty fine ranks participate with zero
counts, and file size is set exactly for shrink/overwrite. Root no longer gathers
whole levels for checkpoint writing. Checkpoint reading remains root-based.

Before opening output, ranks agree on path, complete metadata and ownership.
The writer checks all local values, positive temperatures, the existing EOS
consistency rule, and global restriction error using the original
`64*epsilon` threshold. Numerical state is read-only. Invalid numerical or
metadata input leaves existing files untouched. I/O failure after opening may
leave a partial file, as with the previous writer; no file transaction or
production scaling claim is made.

## Checks

Reused GNU Fortran 13.3 / OpenMPI 4.1.6 Debug build
`build/spray-mpi-migration-gnu-debug`, with bounds/FPE checks. Added a source to
the existing MPI library; compiler, MPI and backend options were unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_amr_reactive_3d pelef_mpi_amr_reactive_3d test_mpi_amr_regrid_3d test_amr_reactive_3d_checkpoint --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(mpi_amr_regrid_3d_np(4|8)|unit_amr_reactive_3d_checkpoint)$' --output-on-failure
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^regression_mpi_amr_balanced_boundary_touching_(dynamic_)?restart_3d$' --output-on-failure
python3 build/amr-mpi-checkpoint-output-checks/compare_parent.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The 4/8-rank unit runs compare entire checkpoint files against serial output
through ratios 2/3, unequal partitions, empty fine ranks, patch movement and
shrink/overwrite. Inconsistent temperature, coarse/fine synchronization and
rank-local time metadata reject without altering the existing checkpoint.
The existing serial checkpoint unit covers accepted and rejected checkpoint
contexts. The application regressions cover static/dynamic reacting physical
patches, load balancing and 4-rank stop / 2-rank restart, with exact serial CSV
comparisons and unchanged conservation gates.

The saved parent serial executable and the refactored serial writer produce
identical periodic checkpoint and coarse/fine CSV bytes for the one-step
entropy-wave checkpoint fixture. This independently checks metadata extraction.
An initial 4/8-rank unit attempt exposed a trailing separator in MPI cell rows
(2 failures); row width was corrected and the focused checks rerun. The earlier
application pass is not counted as final-source verification.

Logs, parent executable, comparison script and source hashes:
`build/amr-mpi-checkpoint-output-checks/`.
Restart artifacts:

- `build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_restart_3d/run-g294jdhz/`
- `build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_dynamic_restart_3d/run-x_h1lnds/`

## Source SHA-256

```text
1cc6cf44b327b267a481a450b9071068177828a250f9c0df4fc46a9fb701e05d  src/parallel/mpi_amr_checkpoint_3d_mod.F90
b65101274f1b111a1647206370abd223742bc28f0dff7d2b5fbd2f6d60dfb7c6  src/io/amr_reactive_3d_checkpoint_mod.F90
85f56ec634f16f7ee0e1ed6dfed4f15df4b84c376d07d25f0c685f809fd6a39a  src/driver/mpi_amr_reactive_3d_application_mod.F90
45dbf1d84c7892a3532d66d2e013c6b3895192e230695a245ee9d4580de8238b  tests/mpi/test_mpi_amr_regrid_3d.F90
```

Final focused checks: **5/5 CTest tests passed, 0 failures, 0 skips**
(`unit-final.log`, `ctest-fixed.log`). The parent-binary byte comparison also
passed (`parent-fixed.log`). Build, provenance, project contract, Python AST
syntax for four verification scripts, and whitespace checks passed.
Subsequent application indentation and an MPI comment changed no executable
statements; the final unit run also asserts the synchronization-specific
rejection diagnostic. No full configuration matrix, installed-binary audit,
external solver parity or physical validation was run.
