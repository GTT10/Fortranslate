# Coupled MPI EB spray application restart — 2026-09-18

Implementation `fcdaf15b6e09e820ce2e7e9741e7b3fa9e506f5a`, base `c4aa13f`.

Fixed and selected MPI EB spray applications now connect the gas snapshot and
particle payload with a shared generation token. A step-suffixed prefix is
reserved with an exclusive `.pending` file. Existing manifest/gas/particle paths
are rejected before writing. The final manifest is linked only after both
payloads close successfully. The original normalized spray namelist and nodal
values bind the particle context; existing gas configuration, mechanism,
geometry and policy checks remain in force.

The gas archive stores a snapshot with its own gas inventory reference and zero
exchanges. The prefix separately preserves the true initial gas inventory and
all accumulated coupled boundary, wall, chemical and injected exchanges. This
keeps the existing gas archive conservation checks meaningful without treating
particle feedback as gas-only history. Restore stages both gas and particles,
checks generation/context, time, step, event, watermark, injected ledger and
combined budget, then publishes all state together. Checkpoints precede events
due at their accepted step's new time; continuation processes pending events
once, including ties. Output, final time and stop/save controls are not included
in the spray context. Regridding/rebalancing remain explicitly unsupported.

## Independently completed checks

GNU Fortran 13.3.0 `/usr/bin/gfortran`, Open MPI 4.1.6 `/usr/bin/mpifort` and
`/usr/bin/mpiexec`, Debug/Ninja. Existing fixed tests-enabled and selected H2/O2
tests-disabled build configurations were reused unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_restart$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-restart-checks/selected-final --selected --scheduled --restart
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, test **48.85 s**, total **48.89 s**. Selected runtime
script passed. Provenance, syntax, contract and whitespace checks passed. No
failures or skips. Initial fixed/selected runs passed before the existing-payload
reservation check and overwrite regression were added; these are revisions of
the same test, not extra distinct evidence.

The scheduled 20 ns fixture saves on two ranks at step two with four live
particles, before the next event. Four-rank continuation produces exactly equal
parsed gas and particle CSV records to uninterrupted four-rank execution. All
initial/final and cumulative ledger columns agree at relative 3e-12, floor 1e-30.
The existing scheduled checks cover tied, zero, final-time and future events.
The test rejects overwrite and verifies all three published files remain
byte-identical. Changed liquid injection temperature rejects the context, and
truncated particle payload rejects before output. The valid payload is restored
after the corruption test. Fixed gas uses elementary H2/O2 initialization;
selected gas is pure N2. Implicit chemistry and all molecular transport are on.

Artifacts: `build/spray-mpi-eb-restart-checks/{build-final.log,build-selected-final.log,ctest-final.log,selected-final.log}`.
Selected final cases: `selected-final/eb_spray_xlg3hftu/` beneath that root.
CTest reports its retained temporary fixture directory.

The application budget retains its existing scale floor of one; this test is
continuation equivalence, not a new tight small-mass conservation qualification.
Manifest publication does not promise fsync/power-loss durability. Context is
bounded to 130000 characters, plus generation; larger checkpoint contexts reject.
No full matrix, standalone gas-only regression, installed-binary audit, scaling,
external solver parity or physical validation was run in this increment.
