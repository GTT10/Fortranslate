# Fixed-mechanism reacting AMR restart — 2026-09-17

Source: parent `cdbc9ebad498cbb5e2027a3b06db2f3e8a2df5b4` plus this change.
Logs and source hashes: `build/amr-fixed-chemistry-restart-checks/`.
Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build, with
`/usr/bin/mpiexec.openmpi` as launcher.

Fixed elementary/full-H2/O2 chemistry no longer disables AMR checkpoint/restart.
The checkpoint layer reconstructs the canonical compiled fixed mechanism and
serializes all reaction records, initial mole fractions, chemistry tolerances and
the actual default integrator policy. Full H2/O2 identity also includes its
generated source SHA. It checks that context before committing any restored
state. This does not permit user-supplied reactions through the fixed path.
Selected contexts keep their existing markers, bundle SHA validation and layout.
Species thermodynamics, transport, geometry and boundary ledgers retain their
existing transactional validation.

New fixed-reacting schemas add 16 to the corresponding fixed nonreacting schema:
18/20 periodic static, 22/24 periodic dynamic, 26/28 physical static and 30/32
physical dynamic (hydro/transport in each pair). Existing schemas are unchanged.
Both serial and MPI applications use the shared checkpoint implementation.

The focused application regression exercises full H2/O2 implicit chemistry with
transport and dynamic physical-boundary AMR, plus elementary explicit chemistry
without transport and with a stationary wall. It compares a serial reference
against a four-rank stop and two-rank restart, byte-for-byte on both level CSVs.
An inactive-chemistry control must differ from each reacting reference. It also
rejects changed chemistry tolerances, integrator records, reaction rates and
model/source identities before output. Boundary-ledger comparisons retain the
existing new-MPI-test tolerance. The unit fixture checks periodic schema 18,
roundtrip and complete state rollback after an integrator mismatch; existing
selected/nonreacting checkpoint tests remain in the same unit executable.

The initial application attempt completed full-H2/O2 dynamic serial/MPI parity,
then rejected the elementary fixture because transport was disabled while wall
motion remained enabled. Its log is retained as `ctest-fixture-failure.log`.
The elementary fixture now uses a stationary adiabatic wall. The reacting
regression uses final time 1e-10 s and CFL 0.0002 to retain two accepted coarse
steps while reducing chemistry work; solver tolerances are unchanged.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_amr_reactive_3d pelef_mpi_amr_reactive_3d test_amr_reactive_3d_config test_amr_reactive_3d_checkpoint --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(regression_mpi_amr_fixed_chemistry_restart_3d|unit_amr_reactive_3d_(config|checkpoint))$' --output-on-failure
python3 -m py_compile tools/check_amr_physical_boundary_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 3/3 tests passed with no skips (application 68.98 s, config 0.00 s,
checkpoint 0.10 s; total 69.11 s). Provenance, Python syntax, project contract
and diff checks passed. Both reacting references completed two coarse steps.
Final run artifacts:
`build/spray-mpi-migration-gnu-debug/tests/mpi_amr_fixed_chemistry_restart_3d/run-e61ev99u/`.

This is bounded numerical restart evidence. It does not establish general
multilevel/multipatch AMR, physical ignition validation or scaling performance.
