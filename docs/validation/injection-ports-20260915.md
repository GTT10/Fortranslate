# Multiple injection ports — 2026-09-15

Source: parent `0a0ee0d` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`.

## Implemented scope

The existing deterministic cone and mass-CDF generator now accepts up to
16 independently positioned/oriented ports sharing a schedule and material
state. Positive mass fractions must close to one; roundoff is normalized.
The total initial mass and mass flow are distributed among ports. Each port
creates `parcel_count` parcels, with ordered consecutive IDs. A failed later
port leaves no partially generated output; the application appends only a
complete event. Storage checks and restart parent/event counts use the full
number of parcels per event.

Active geometry and mass shares participate in checkpoint identity. Count
zero preserves the previous single-port path/context. In multiple-port mode,
unused single-port coordinates are not required to lie inside the domain.
The core generator validates all active axes; application validation also
checks positions and rejects nonzero unused array entries.

Ports share temperature, speed, cone angle, size distribution, liquid and
schedule. This does not add independent injector schedules/fuels, collisions,
resolved nozzle flow or MPI/AMR particle support.

## Focused evidence

- Three core cases passed: ports, injection and size distribution. The new
  case checks two opposing ports with 25/75 mass shares, six ordered IDs,
  correct weighted total momentum, per-port mass/position, fraction closure
  and failure of a later port without partial output.
- Thirteen CLI tests passed in the initial 10.16 s run. The new test checks
  initial plus scheduled injection (16 parcels), a 3.5e-11 kg injected-mass
  ledger, conservation, exact gas/parcel restart, changed share rejection,
  invalid fractions and combined storage limit rejection. The CLI suite was
  rerun and passed in 10.23 s after making the single-port domain check conditional; its multiple-port
  case now uses a 1 mm domain with valid active ports and an unused default
  single-port center outside that domain.
- The supplied 4^3 example completed 20 steps to 2e-6 s, with four scheduled
  events plus initial injection across two ports, ending with 40 live parcels.
  Injected mass was 1.0000000000000007e-10 kg. Maximum relative system mass and
  energy errors were 6.4208978691e-16 and 1.6466275684e-15.
- Provenance preflight, project contract, Python syntax and whitespace passed.
  This is integration/conservation evidence, not a physical injector experiment.

Artifacts: `build/water-ports/` contains the input example's outputs, checkpoint,
log and `result.json`; build logs are `build/port-spray-debug/ports-build.log`
and `ports-final-build.log`. CLI cases use temporary directories.

```sh
cmake --build build/port-spray-debug --target test_spray_core pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^spray_(ports|injection|size_distribution|application_cli)$'
mkdir -p build/water-ports
cd build/water-ports
../port-spray-debug/pelef_spray_3d ../../cases/spray_3d/water_ports.nml
```

Use a fresh output directory. Additional checks from the repository root:

```sh
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
