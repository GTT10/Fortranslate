# MPI boundary-touching AMR application/restart — 2026-09-17

Source: parent `0e19d70` plus this change.

The fixed and selected MPI frontends now accept static patches touching
nonperiodic physical faces. The application uses valid geometry rather than
requiring an interior patch. Physical-boundary regridding remains explicitly
rejected before initialization. Interior dynamic cases retain their existing
behavior.

For isothermal walls, each rank computes the coarse wall constraint and, when
it owns fine cells on a physical wall, the fine constraint at refined spacing.
The fine bound is converted through the ratio-squared transport subcycle count
before the global minimum is taken. Ranks without fine cells participate in
the reduction. Boundary transfer and existing geometry/chemistry-aware
checkpoint records use the already connected shared implementation; no schema
layout changed.

The focused regression runs full-H2/O2 reacting moving-wall/transport and
elementary reacting corner cases. It compares serial uninterrupted execution
with four-rank checkpoint creation and two-rank restart, including coarse/fine
CSV bytes and boundary ledgers at the established MPI tolerance. Existing
boundary/chemistry context mismatch probes and inactive-chemistry comparison
remain included. The selected frontend is wired but was not independently
built or qualified in this check.

Commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_amr_reactive_3d pelef_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^regression_mpi_amr_boundary_touching_restart_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build. Logs and source hashes:
`build/amr-mpi-boundary-touching-app-checks/`.
Result: 1/1 test passed in 54.40 s (CTest total 54.43 s), no failures or
skips. Provenance, project contract and diff checks passed. No Python source
changed. Final application artifacts:
`build/spray-mpi-migration-gnu-debug/tests/mpi_amr_boundary_touching_restart_3d/run-6i5g_s64/`.

This is application/restart evidence, not external solver parity, production
scaling or physical validation. MPI boundary-touching dynamic remap and general
multipatch/multilevel hierarchies remain open.
