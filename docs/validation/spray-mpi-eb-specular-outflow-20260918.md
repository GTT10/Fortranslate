# MPI specular EB source outflow — 2026-09-18

Implementation `12b9884`, base `93f1736`.

The MPI source `eb_outflow_change(nvar)` now also works without impact films.
The specular drift evaluates mixture liquid energy, clips EB/exterior flight,
zeros escaped component inventory and reduces local exterior losses once across
ranks. Accepted steps accumulate global component-to-species inventories while
keeping EB wall exchange separate. Escaped records are compacted after each
accepted substep before any subsequent gas sampling. Final ownership migration
and whole-interval rollback remain inside the transaction.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Selected library compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_step_(2|4)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**, 0 skips, 2.10 s total, configured `/usr/bin/mpiexec`.
Builds, provenance, syntax, contract and diff checks passed. Initial tests failed
2/2 at the new reflection/escape assertion: retaining escaped records until final
migration caused later source work to fail. Moving compaction into accepted
substeps resolved the failure. Initial logs are retained; no tolerance changed.

A single 0.25/0.75 mixed parcel reflects from the embedded wall and exits the
root during a 1 ms interval, with drag/heat/evaporation/breakup disabled. Final
population is empty and watermark stays 1. Wall y momentum is positive, exterior
mass negative, and mapped species follow mixture fractions to 1e-23 kg. Global
composite inventory closes against wall plus exterior at existing relative
3e-12 tolerance with 1e-30 floor. A one-substep cap preserves original parcel/
component payload, watermark and both ledger sentinels. Existing impact/film,
KH, transfer and ownership checks remain in the executable.

Artifacts: `build/spray-mpi-eb-specular-outflow-checks/` contains `provenance.log`,
`build.log`, `ctest.log`, `build-selected.log`, `build-final.log`,
`ctest-final.log`, and `build-selected-final.log`.

Specular outflow is not yet exposed by the MPI reactive wrapper/application.
No physical validation or external-solver parity is claimed.
