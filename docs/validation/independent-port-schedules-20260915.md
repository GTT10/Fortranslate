# Independent injection-port schedules — 2026-09-15

Source: parent `1fc5220` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`.

## Implemented scope

Each configured port can now have its own mass flow, start/end time and event
interval. The scheduler merges the individual event streams by timestamp and
port index, handles simultaneous events before the flow step, and limits the
step to the next event. Initial mass allocation remains separate and still
uses port fractions; scheduled rates are per-port absolute kg/s values.
Per-event parcel count remains shared. Geometry, fuel/temperature/size/speed
behavior uses the existing injection generator.

The checkpoint's aggregate event count is sufficient to replay the deterministic
merged schedule. Restart rebuilds per-port counts and checks time consistency,
injected mass and parent-parcel inventory. A schedule change rejects restart
through context identity. The existing shared scheduler remains the default.
Active independent schedules require zero global mass flow; unused array data,
invalid intervals/windows and more than one million total events are rejected.
Replay costs O(events*ports); scalability of million-event restarts is not claimed.

## Focused evidence

- All fourteen CLI tests passed in 10.93 s. The new case uses two ports starting
  together at 1e-7 s, with intervals 2e-7/3e-7 s, ends 5e-7/6e-7 s and rates
  1e-4/2e-4 kg/s. Four individual port events create eight parcels and inject
  1.4e-10 kg. Exact IDs and velocity signs verify timestamp/port ordering.
- Stops after steps 1 and 3 exercise a pending simultaneous event and a later
  pending event. Both resumes reproduce uninterrupted gas and parcel files
  exactly. The suite also checks coupled balances, changed rate rejection,
  zero active interval rejection and incompatible global/independent flow.
- The supplied 4^3-cell example completed 20 steps to 2e-6 s. Eight additional
  port events plus initial injection end with 40 parcels and an injected-mass
  ledger of 3.5999999999999984e-10 kg. Maximum relative system mass and energy
  errors were 5.0413044091e-16 and 3.3415787924e-15.
- Provenance preflight, project contract, Python syntax and whitespace passed.
  No new mechanism, MPI/AMR particles or physical injector experiment is claimed.

Artifacts: `build/water-independent-ports/` contains the example's outputs,
checkpoint, log and `result.json`; build log is
`build/port-spray-debug/port-schedules-build.log`. CLI cases use temporary directories.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^spray_application_cli$'
mkdir -p build/water-independent-ports
cd build/water-independent-ports
../port-spray-debug/pelef_spray_3d ../../cases/spray_3d/water_independent_ports.nml
```

Use a fresh output directory. Additional checks from the repository root:

```sh
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
