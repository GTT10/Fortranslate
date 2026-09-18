# MPI EB impact interval with exterior outflow — 2026-09-18

Implementation `5e48526`, base `113cb49`.

`advance_mpi_level_set_impacts` accepts optional `outflow_change(5+nc)`.
Ranks collectively check presence and output size. Each local trajectory is
clipped to its first exterior crossing before the collective EB contact batch.
An earlier EB contact's fractional time is converted back to the full remaining
flight interval. An exterior-first particle contributes its represented mass,
momentum, liquid internal plus kinetic energy and component masses to a private
negative system-change ledger, then becomes inactive. Final compaction removes
it without recycling IDs. MPI_SUM gives every rank the same exterior inventory.

Exterior events count against the existing per-lineage hit limit. A failed
rank agreement, impact or event limit publishes neither populations nor ledgers.
Wall receipt remains separate. The tree migration wrapper and MPI source,
reactive and application paths are not yet connected to this option.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja, reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Selected library compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_(2|4)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**, 0 failures/skips, 0.98 s total, using the configured
`/usr/bin/mpiexec`. The added fixture places one water parcel only on the last
rank. One splash produces four children; all escape within 0.1 s. Global hits
are five and the watermark advances by four. Film plus wall receipt minus
initial inventory closes against exterior inventory with the existing relative
3e-12 tolerance and 1e-20 floor. Component loss equals mass loss.

A hit limit of one restores the pre-flight position, films, watermark and both
ledger sentinels. Supplying the optional output on only some ranks rejects
collectively. Existing contact, repeated-rebound and tree-migration checks remain.
An initial run before adding the mismatch check also passed 2/2; it is not an
additional distinct test count. No external-solver or physical validation ran.

Artifacts: `build/spray-mpi-eb-outflow-checks/` contains `provenance.log`,
`build.log`, `build-final.log`, `ctest.log`, `ctest-final.log`, and
`build-selected.log`. Builds retain existing compiler warnings.
