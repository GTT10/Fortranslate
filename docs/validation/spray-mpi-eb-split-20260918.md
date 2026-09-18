# MPI EB computational splitting — 2026-09-18

Implementation `0198749`, base `3d7b023`.

The reactive wrapper connects the existing collective splitter after both
particle halves and before publication. The existing maximum-multiplicity option
controls it; no new input control is introduced. Capacity failure rolls back gas,
temperature, particles, films, shared IDs and budgets. The application records
cumulative split IDs in the gas checkpoint prefix and validates the remaining
watermark against injection and four-child KH/splash groups on restore.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja; reused fixed tests-enabled and selected H2O2 tests-disabled builds.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_reactive pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_reactive_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_impact.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-split-checks/runs --split
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final reactive CTest **2/2 passed**, total **4.83 s**. The fixture checks late
capacity rejection and full rollback, then five descendants with multiplicity
20, retained owner/composition/lineage, and composite conservation within relative
3e-12, scale floor 1e-30. An initial fixture failed both ranks' tests because
disabling KH children left a nonzero KH reference parameter. The corrected ETAB
fixture resets that parameter; no production tolerance was changed.

The application driver passed five launches: full four-rank, stopped two-rank,
resumed four-rank, rejected context mismatch, and rejected altered split count.
Two initial water parcels split into four records with IDs 1–4 and multiplicity
below 60. The step-two V1 archive carries watermark 4 and cumulative split count
2. Changing that count to 1 rejects restore before output. Gas, particles and
budgets match changed-rank continuation within relative 3e-12, floor 1e-30.
Prescribed regrids remain active; accepted rebalance movement is not asserted.
New final-time injections wait until a subsequent reactive update to split.

Both executable builds, provenance, Python syntax, project contract and whitespace
checks passed. Selected execution, broad matrix, installed audit, external solver
parity and physical validation were not run. No final checks were skipped.
The application test is registered as `spray_mpi_eb_split_application` and was
invoked directly for this record.

Artifacts: `build/spray-mpi-eb-split-checks/`, including `provenance.log`,
`build-final.log`, `ctest-final.log`, `application.log`,
`build-selected-final.log`; initial failed CTest in `ctest.log`.
Inputs, outputs and checkpoints: `runs/impact_8_sl7bz7` below that directory.
