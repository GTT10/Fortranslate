# MPI EB acceleration application and restart — 2026-09-18

Implementation `c1a13e9`, base `e7a6c1b`.

The MPI EB spray namelist exposes bounded `acceleration(3)` and forwards it to
both reactive wall-model branches. Cumulative force is separate from the common
boundary/wall/chemical record, enters composite and restored budget checks, and
appears as a `force` column in MPI spray budget CSVs (zero when disabled).
Nonzero acceleration stores an additional force row in the coupled gas prefix;
restart context binds the input/mode. Reading stages and checks finite values,
then publishes force only after gas and particle/film restoration succeeds.
The particle format is unchanged. Zero-force checkpoints omit the extra row.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Selected executable compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_force_application$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, 0 skips, 19.52 s total with configured
`/usr/bin/mpiexec`. Four launches cover 4-rank full, 2-rank checkpoint, 4-rank
resume and changed-context rejection. The deposition/film case uses acceleration
(0,0,1e6) m/s², three regrids and rebalancing. Positive z force is checked.
Gas, particle, film and budget rows, including force, match full/resumed output
at existing relative 3e-12 comparison and 1e-30 floor.

Initial fixed and selected builds failed on three Fortran lines exceeding 132
columns. Wrapping those lines fixed both builds. A test launched after the failed
fixed build used the stale executable and failed; it is excluded from current
implementation evidence. Initial logs remain. Final builds, syntax, provenance,
contract and diff checks passed. No tolerance changes were made.

Artifacts: `build/spray-mpi-eb-force-app-checks/` contains `provenance.log`,
initial `build.log`, `build-selected.log`, `ctest.log`, and final
`build-final.log`, `build-selected-final.log`, `ctest-final.log`.
Final fixture: `build/spray-amr-reactive-debug/tests/spray_mpi_eb_force_application/impact_zlivrrwp/`.
This is application/restart evidence; broader physical and external-solver
qualification remains open.
