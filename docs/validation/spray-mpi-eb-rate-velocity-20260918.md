# MPI EB profile-derived injection speed — 2026-09-18

Implementation `ec44b31`, base `363be86`.

Ports can derive injection speed from their scaled profile integral divided by
event duration, mixture density, outlet area and discharge coefficient. This
uses the same mass-flow relation as the regular spray application's rate-velocity
path. Cone/swirl directions are retained. Both particle payloads and injected
momentum/energy use the computed velocity. Zero-mass events reserve no IDs.
Settings and profile intervals remain in normalized checkpoint context.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja; reused fixed tests-enabled and selected H2O2 tests-disabled builds.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_nozzle.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-rate-velocity-checks/runs --rate
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The driver passed five MPI launches: full two/four-rank runs, two-rank stop,
four-rank resume, and intentional negative-radius rejection. No unexpected
failures or skips. Both executable builds, provenance, Python syntax, contract
and whitespace passed. Selected execution was not run. CTest registration is
`spray_mpi_eb_rate_velocity_application`; this execution called the driver directly.

Two finite solid/hollow ports retain swirl and mass/number size weighting. They
use constant and linearly increasing profiles with coefficients 0.8 and 0.5.
A zero-scale event at 10 ns consumes no IDs. Final-time event doses are 5e-12 and
2.125e-11 kg; measured fresh-particle speeds match the independent formula within
relative 1e-12. Event masses match within relative 1e-13. Position, swirl,
size bounds, weighting and sequential IDs retain the prior nozzle checks.
The application's coupled budget checks include injected momentum and energy.
Gas and particle CSVs agree across two/four ranks and after restart within
relative 3e-12, scale floor 1e-30. Negative radius rejects before output.

Artifacts: `build/spray-mpi-eb-rate-velocity-checks/{build.log,build-selected.log,application.log,provenance.log,contract.log}`;
per-run inputs, outputs and checkpoints under `runs/nozzle_2yvfzuvg` there.
No broad matrix, installed audit, external solver parity, physical nozzle-flow
validation or scaling check was run. The profile's entire interval mass is still
injected as a packet at its start; this is not continuous-time inlet integration.
