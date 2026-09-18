# MPI AMR spray force application and restart — 2026-09-18

Implementation `18c85668a4a0cf6af88fe34bf8cd08497b4bbf7b`, based on `cba4fad`.

Fixed/selected MPI AMR spray applications read optional `&amr_spray_force`
acceleration, validate its finite component bounds and bind nonzero values to
restart identity. Accepted coupled steps accumulate a separate source ledger,
reported as `Source inventory:`. Both film and ordinary branches pass the force
controls to the coupled API.

`mpi_amr_spray_checkpoint` accepts optional `source(nvar)`. Presence agrees
collectively; writes require finite rank-identical values before I/O. A gas
context marker binds the optional `SPRAY_SOURCE_V1` particle-payload extension.
Reads validate the extension before migrating candidate parcels and publishing
gas, particles and source together. Default-disabled application runs omit the
extension and retain the previous format. Immutable publication is unchanged.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `pelef_mpi_amr_spray_3d` and selected H2O2
`pelef_mpi_amr_spray_3d_selected` built; selected executable was not run.
CTest **2/2 passed**, 9.55 seconds, zero failures/skips:
`spray_mpi_amr_force_application` (4.49 s), `spray_mpi_amr_restart` (5.02 s).

The force fixture uses acceleration `[100,-200,300]`, scheduled coarsening and
injection. Continuous 2-rank output matches 2-to-4 rank checkpoint continuation
for gas, parcels and cumulative force inventory within the existing `2e-12`
relative / `1e-25` absolute tolerances. Momentum-source signs are checked and
the saved payload contains the source record. Changing acceleration rejects
restart before CSV publication. Existing context/truncation/overwrite checks
also pass. The default-disabled restart fixture passes separately. An indentation
correction after these checks changed no Fortran tokens.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests or tolerance changes. Force with films, LES, physical
boundaries or chemistry was not separately run at application level. No full
matrix, external solver parity, scaling or physical validation was performed.

## Commands and artifacts

Logs: `build/spray-amr-force-application-checks/{provenance,build,build-selected,ctest}.log`.
Run inputs, logs and CSV/checkpoint artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_{force_application,restart}/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(force_application|restart)$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
