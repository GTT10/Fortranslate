# MPI AMR cell-count balancing — 2026-09-17

Source: parent `a416453` plus this change. Logs and source hashes:
`build/amr-mpi-load-balance-checks/`.
Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build with bounds checks.

The optional `load_balance = .true.` AMR input invokes a deterministic contiguous
slab partitioner after initial/restart distribution and after geometry changes.
Coarse cells and fine children each count once. A capacity search minimizes the
maximum rank cell count subject to at least one coarse column per rank. An
equally good current layout is kept; only strict improvement invokes the direct
MPI repartition API. Rank policy is checked collectively before branching.
The default remains disabled. The shared serial reader accepts the option but
serial execution needs no repartitioning.

The option is independent of checkpoint physics and is not serialized in the
physical/dynamic context signature. Restart reconstructs the global hierarchy,
then chooses ownership for the current rank count. No checkpoint schema changed.
This is a cell-count proxy, not measured chemistry/transport work: stiffness,
subcycling and per-cell runtime are not modeled.

The four/eight-rank remap unit checks the partitioner's maximum cell count
against exhaustive enumeration of all contiguous partitions of its small grid,
exact preservation of serial states and temperatures, and no repeated transfer
for an already optimal layout. Existing explicit repartition/remap and rollback
checks remain. Two application regressions use an eight-column coarse grid,
reacting moving-wall/transport and elementary corner cases, static and dynamic
patches, and four-to-two-rank restart. They require actual balancing, exact
serial coarse/fine CSVs and the established boundary-ledger tolerance. Dynamic
cases require another balance after patch expansion. Existing context rejection
and inactive-chemistry probes remain included.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_regrid_3d pelef_mpi_amr_reactive_3d pelef_amr_reactive_3d --parallel 4
python3 -m py_compile tools/check_amr_physical_boundary_restart.py
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(mpi_amr_regrid_3d_np(4|8)|regression_mpi_amr_balanced_boundary_touching_(dynamic_)?restart_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 4/4 tests passed, no failures or skips (unit four ranks 1.37 s,
eight ranks 1.10 s; balanced static restart 54.42 s; balanced dynamic restart
83.28 s; total 140.21 s). Provenance, Python syntax, project contract and diff
checks passed. The build succeeded with real-equality warnings.

The static wall and corner fixtures chose coarse widths `1 1 2 4` and
`1 1 1 5`, respectively. Dynamic expansion changed both to `2 2 2 2`.
Final application artifacts:
`build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_restart_3d/run-nbk9ybz9/`;
`build/spray-mpi-migration-gnu-debug/tests/mpi_amr_balanced_boundary_touching_dynamic_restart_3d/run-j26edrha/`.

General patches/levels, distributed I/O, measured performance/scaling and
physical qualification remain open. The selected frontend shares this path
but was not independently built or qualified here. No full configuration
matrix or external solver parity run was performed.
