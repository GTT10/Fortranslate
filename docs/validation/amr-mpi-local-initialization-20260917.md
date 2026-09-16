# Rank-local MPI AMR initialization — 2026-09-17

Parent: `6b0364a` on local `main`.

Fresh two-level MPI AMR runs now evaluate the shared analytic problem initializer
on each rank's coarse/fine coordinates. Domain extents stay global, preserving
entropy-wave phase and periodic hotspot distances. Only local array dimensions
change. Fine ranks with zero cells skip evaluation. Parent cells are averaged
in the same child order as the serial implementation; coarse temperatures are
then recovered locally. Initial conserved and elemental inventories use the
distributed diagnostic reductions. Reduction roundoff can differ from serial;
conservation thresholds are unchanged.

Fresh runs no longer allocate full state/temperature arrays on root or scatter
initial fields. Coordinate vectors remain replicated. Checkpoint reads and writes
still gather/scatter complete fields and retain their existing schemas. This is
analytic initialization, not a distributed arbitrary-file reader or a production
scaling measurement. No broader AMR or physical qualification is implied.

## Reproduction and artifacts

Reused `build/spray-mpi-migration-gnu-debug`, GNU Fortran 13.3 / OpenMPI 4.1.6,
Debug with bounds/FPE checks. No compiler/MPI/backend configuration changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_amr_reactive_3d pelef_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^regression_mpi_amr_balanced_boundary_touching_(dynamic_)?restart_3d$' --output-on-failure
python3 build/amr-mpi-local-initialization-checks/check_initial.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The local check script derives inputs from `amr_entropy_wave_plm.nml`, sets
`final_time=1e-12`, refinement ratio 3, and each of `entropy_wave`,
`uniform_reactor`, `reactive_hotspot`. It runs the serial and 8-rank applications
and compares both CSV files byte for byte. All three comparisons passed,
including one coarse column per rank and empty fine ranks. The first script
invocation incorrectly passed an output prefix to the serial CLI and exited
with its usage error; outputs were moved into namelist settings and the complete
three-case check was rerun successfully.

Logs, generated inputs, CSVs and script: `build/amr-mpi-local-initialization-checks/`.
Restart artifacts:

- `build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_restart_3d/run-alcnx5ig/`
- `build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_dynamic_restart_3d/run-mbdkd0rw/`

Source SHA-256 (`src/driver/mpi_amr_reactive_3d_application_mod.F90`):
`aab4510a43191639bd840aefac4eae45a11808bac3a4f8afe6019974b73a3f70`.

Completed CTest regressions: **2/2 passed, 0 failures, 0 skips** (static and
dynamic balanced physical-boundary restart). Both compare serial uninterrupted
CSV with 4-rank stop / 2-rank restart, for reacting wall-transport and elementary
corner cases. The three additional initialization comparisons also passed.
Build, provenance preflight, project contract, Python AST syntax (four check
scripts) and diff whitespace checks passed. No full configuration matrix,
installed-binary audit, external solver parity or physical validation was run.
