# Per-port injection conditions — 2026-09-15

Source: parent `aa3438b` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`.

## Implemented scope

Active ports may override speed, diameter scale, temperature and cone angle.
The -1 sentinel inherits the corresponding global input. Unused values must
remain -1; the resolved state is validated through the existing injection
kernel. Initial, shared-schedule and independent-schedule injection use these
conditions. Different single-drop masses are represented by different
multiplicities while preserving each port's requested mass allocation.

Overrides are part of checkpoint identity; no-override contexts remain.
The liquid and diameter-distribution shape/truncation bounds are still shared.
This does not add different fuels per port, mixture droplets or physical
nozzle-flow qualification.

## Focused evidence

- Three core cases passed: port states, ports and injection. The new case
  checks 10/20 m/s opposing ports, 10/20 micrometer diameters, 300/350 K
  temperatures, their total thermal/kinetic energy and invalid override dimensions.
- Fifteen CLI tests passed in 11.65 s. A new shared-schedule case verifies
  all initial and later parcel speeds, temperatures and diameters, conservation,
  exact restart, changed-speed rejection and an invalid later-port diameter.
  The independent-schedule test now uses different speeds and checks exact
  velocity/ID order through both restart points. Four CTest entries passed
  in 11.69 s.
- The supplied combined case completed 20 steps to 2e-6 s with 40 parcels.
  It uses independently scheduled ports at 80/120 m/s, 20/40 micrometers,
  310/340 K and 0.05/0.15 rad cones. Injected mass was
  3.5999999999999984e-10 kg. Maximum relative system mass and energy errors
  were 7.2520837694e-16 and 4.1642485107e-15.
- Mechanism provenance, project contract, Python syntax and whitespace passed.
  These are input/integration/conservation checks, not physical spray validation.

Artifacts: `build/water-port-states/` contains the example output, checkpoint,
log and `result.json`; build log is `build/port-spray-debug/port-states-build.log`.
CLI fixtures use temporary directories.

```sh
cmake --build build/port-spray-debug --target test_spray_core pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^spray_(port_states|ports|injection|application_cli)$'
mkdir -p build/water-port-states
cd build/water-port-states
../port-spray-debug/pelef_spray_3d ../../cases/spray_3d/water_port_states.nml
```

Use a fresh output directory. Additional checks from the repository root:

```sh
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
