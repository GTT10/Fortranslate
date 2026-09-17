# MPI EB spray application regrids and continuation — 2026-09-18

Implementation `d3e573ad028a4d1a7e09bbbe97d1cf9ff7daf4b2`, base `295fb96`.

The fixed/selected MPI EB spray applications now call the coupled EB remap
when activating a regrid stage. Gas, geometry, topology, owners, particles and
component masses change together. Existing prescribed schedules and the existing
temperature-adaptive tagger feed the same activation path. The original nodal
field is retained at the maximum schedule lattice resolution.

On spray checkpoint restore, the policy reader records the saved stage without
activating it globally. A private hierarchy/partition/geometry is prepared for
that stage using the existing gas remapper. Gas and particle archives are read
into private state and their generation, metadata, stage/time and combined
inventory are checked before publishing the new hierarchy and state. Ordinary
gas-only policy readers retain their previous activation path.

## Completed focused checks

GNU Fortran 13.3.0 `/usr/bin/gfortran`, Open MPI 4.1.6 `/usr/bin/mpifort` and
`/usr/bin/mpiexec`, Debug/Ninja. Existing fixed tests-enabled and selected H2/O2
tests-disabled configurations were reused unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_regrid_application$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-regrid-app-checks/selected-final --selected --scheduled --restart --regrid
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, test **45.95 s**, total **45.98 s**. Selected application
execution passed. Provenance, Python syntax, contract and whitespace passed.
An initial compile failed because one expanded argument line exceeded Fortran's
line length; it was wrapped without semantic change. Test commands following
that failed build ran old binaries and rejected the newly requested regrid mode;
those results are not counted as validation of this implementation. Final runs
were started only after successful builds. No skips or tolerance changes.

The 20 ns fixture starts with left refinement, moves it right at 5 ns, coarsens
to root-only at 10 ns and refines left again at 15 ns. Scheduled injection has
tied, zero-scale, final-time and pending future events. Both two/four-rank runs
report all three regrid stages and preserve the existing gas/particle rank parity
checks. A two-rank checkpoint at step two stores the root-only stage with four
live particles. Four-rank continuation processes the pending injection and later
refinement, producing exactly equal parsed gas and particle CSV records to an
uninterrupted four-rank run. All initial/final and cumulative ledger columns
agree at relative 3e-12, floor 1e-30. The same script checks immutable publication,
changed spray context and truncated particle rejection.

Artifacts: `build/spray-mpi-eb-regrid-app-checks/{build-final.log,build-selected-final.log,ctest-final.log,selected-final.log}`.
Initial failure logs are retained alongside them. Selected final cases:
`selected-final/eb_spray_m7ystbfa/`; CTest LastTest.log identifies the fixed fixture.

This is prescribed-regrid continuation evidence. Temperature-adaptive tagging
uses the same integration path but was not independently exercised here.
Rebalancing and regridding together were not newly qualified in this fixture.
No full matrix, standalone gas-only regression, convergence, scaling, installed
binary audit, external solver parity or physical validation was run.
