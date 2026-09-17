# Scheduled EB spray injection — 2026-09-18

Implementation: `a5081d44e211938649184927b589dea76e78e5e8`, based on `271f350`.

The shared fixed/selected EB driver accepts up to 1024 nondecreasing event
times and mass multipliers. Each positive event injects the configured parcel
count with scaled represented mass. Event candidates must lie in fluid and fit
the retained population limit before population, component masses, IDs and the
integrated injection ledger commit. Zero-dose events consume no IDs. Timesteps
clip to the next event; ties and events exactly at final time are processed.
Defaults retain the original event at time zero. Nozzle and liquid tables are
shared by all events.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; unchanged existing fixed
tests-enabled and selected H2/O2 tests-disabled configurations.
Compiler `/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(schedule|application)$' --output-on-failure --output-log build/spray-eb-schedule-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/spray-eb-schedule-checks/selected --selected --scheduled
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, total **0.90 s**: schedule **0.53 s**, default application
**0.34 s**. Selected scheduled execution passed independently. Provenance,
syntax, contract and whitespace checks passed. No build/test failures or skips.

Schedule: times `[5,5,5,10,20]` ns, multipliers `[1,0,1,0.5,2]`, two parcels per
positive event. The initially empty population advances gas before injection.
Assertions require eight final parcels with IDs 1–8, identical states for the
two positive tied events, heating/evaporation/reflection in evolved populations,
unevolved final-time parcels, and total injected mass 4.5e-10 kg at relative
1e-14. Every gas record reaches 20 ns. Decreasing event times are rejected, and
the existing input/output alias check still passes. The fixed and selected gas
compositions differ; this is not a fixed/selected parity claim.

Artifacts: `build/spray-eb-schedule-checks/` contains `build-final.log`,
`build-selected.log`, `ctest.log`, `selected.log` and the selected case directory
`selected/eb_spray_u4zcxq2x/`. Fixed CTest cases are under
`build/spray-amr-reactive-debug/tests/spray_eb_schedule/`.

Checkpoint/restart, independent ports, time-varying rate-profile integration,
outer particle boundary events and MPI/AMR EB ownership remain open. No full
configuration matrix, convergence, external parity or physical validation ran.
