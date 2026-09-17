# Scheduled distributed EB spray injection — 2026-09-18

Implementation `a46e4fb519341d2778ff821f82d9147bb74f14b5`, base `f2e6ea3`.

The fixed/selected MPI EB spray applications now accept up to 1024 ordered
injection times and nonnegative mass scales. A replicated cone template fixes
per-drop mass and composition; event scales multiply multiplicity. Only rank
zero appends new records before collective finest-fluid-owner migration.
Candidate populations publish after successful migration. Global live count is
bounded by the configured cap and 125000; IDs increase only for positive events.
The common step clips to the next event. Ties, zero scales, final-time events
and pending events after final time have explicit semantics.

The initial budget is gas before injection, including when the first event is
at zero. An `injected` column records mass, momentum, liquid internal/kinetic
energy and vapor-species inventories. Balance subtracts these source totals.
The gas-only output schema is unchanged. Checkpoint/restart remains rejected.

## Completed focused evidence

GNU Fortran 13.3.0 `/usr/bin/gfortran`, Open MPI 4.1.6 `/usr/bin/mpifort` and
`/usr/bin/mpiexec`, Debug/Ninja. Existing fixed tests-enabled and selected H2/O2
tests-disabled configurations were reused without configuration changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(schedule|application)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-schedule-checks/selected --selected --scheduled
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, **45.94 s** total: default zero-time injection **21.91 s**,
schedule **24.00 s**. Selected scheduled execution passed. Provenance, syntax,
contract and whitespace passed; no failures or skips.

Scheduled tests use times [5,5,5,10,20,30] ns with scales [1,0,1,.5,2,1],
final time 20 ns and two records per positive event. Assertions verify eight
unique consecutive IDs, five processed events, tied-event identical states,
the final two particles at their injection state, a still-pending future event,
and injected mass 4.5e-10 kg within relative 1e-14. Earlier records reflect and
change temperature. Every gas and parcel CSV field agrees between two and four
ranks within relative 3e-12 with 1e-30 floor. Decreasing schedules and unavailable
restart are rejected before output. Fixed gas retains the elementary H2/O2
initialization; selected gas is pure N2. Both enable implicit chemistry and all
molecular transport terms.

Artifacts: `build/spray-mpi-eb-schedule-checks/{build-final.log,build-selected.log,ctest.log,selected.log}`.
Selected output directory: `selected/eb_spray_jk4fjg2a/` under that artifact root.
CTest logs identify retained temporary fixture directories.

This is scheduling and MPI integration evidence, not tight independent
small-mass conservation, physical validation, convergence, scaling or full
configuration coverage. No standalone gas-only regression was repeated for this
increment. Regridding, rebalancing, restart, independent ports and EB liquid
films remain open for this executable.
