# Independent MPI EB injection ports — 2026-09-18

Implementation `646d1ff5d4b9e689d8e864556c4d3fc602b15d48`, base `9a532f1`.

The MPI EB spray application now supports 1--16 cone injection ports. Each
port may specify count, center, axis, speed, half-angle, diameter, temperature,
base mass and composition. Missing entries inherit the existing scalar inputs.
`event_ports` selects a port for each event in the globally ordered table;
equal-time events may select distinct ports. Each port's mixture properties
construct its own injection template. Per-drop component mass persists through
existing coupling; material ID records the injection port. CSV now includes
`material_id` explicitly.

Positive events reserve the selected port's count of IDs and scale its template
multiplicity. Zero events reserve no IDs. Injected mass, momentum, energy and
species use that port's composition and temperature. Restart checks the sum of
counts for processed positive events rather than a single constant event count,
and rejects material IDs beyond the configured port count. Normalized port
settings and event selections participate in the strict checkpoint context.
Existing single-port input still works; pre-extension checkpoint contexts are
not silently upgraded.

## Completed focused checks

GNU Fortran 13.3.0 `/usr/bin/gfortran`, Open MPI 4.1.6 `/usr/bin/mpifort` and
`/usr/bin/mpiexec`, Debug/Ninja. Existing fixed tests-enabled and selected H2/O2
tests-disabled configurations were reused unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_ports|spray_mpi_eb_application)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-ports-checks/selected --selected --scheduled --restart --regrid --ports
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, total **69.00 s**: legacy single-port application **26.16 s**,
new port application **42.80 s**. Selected runtime script passed. Provenance,
Python syntax, project contract and whitespace passed. No failures, skips or
tolerance changes.

The 20 ns prescribed-regrid fixture uses two ports. Port one inherits two 350 K
water records, x=1.5 mm, speed 20 m/s and 10 micrometer diameter. Port two has
three 360 K records, x=3 mm, speed 30 m/s, 20 micrometer diameter, twice the base
mass and a 25/75 water/H2 component split. These are synthetic constant-property
component tables for software integration, not physical liquid hydrogen data.
Event selectors [1,2,2,1,2,1] at [5,5,5,10,20,30] ns use scales [1,0,1,.5,2,1].
Tests require ten consecutive IDs, correct port IDs, component mass sums,
injected mass 7.5e-10 kg (relative 1e-14), reflected earlier particles and the
three final-time port-two records at their untouched injection state/composition.

Two/four-rank gas and particle outputs agree at relative 3e-12 (floor 1e-30).
Two-rank root-only checkpoint to four-rank continuation preserves pending
injection, subsequent refinement and exact parsed gas/particle outputs against
uninterrupted execution. All cumulative ledger columns retain their existing
relative 3e-12 comparison. Existing overwrite, changed-context and truncated
payload rejection checks also run. Fixed gas uses elementary initialization;
selected gas is pure N2, with implicit chemistry and all molecular transport.

Artifacts: `build/spray-mpi-eb-ports-checks/{build-final.log,build-selected.log,ctest.log,selected.log}`.
Selected cases: `selected/eb_spray_jshs7daj/`; CTest LastTest.log identifies fixed
case directories. Independent per-port rate-curve files, EB liquid film impacts,
physical qualification and scaling remain open. No full matrix, installed audit,
external solver parity or physical validation was run.
