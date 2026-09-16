# Shared profile integration and MPI injection — 2026-09-17

Parent `777c1ba` plus this change. Source hashes and MPI test/build/preflight logs:
`build/spray-mpi-profile-checks/`.

A common validated piecewise-linear integrator now serves serial and MPI spray.
MPI profile points cover the complete event window, replacing fixed event mass
with each interval integral. Zero-dose events advance the event counter without
particles or IDs. Profile identity enters the checkpoint context; absent profiles
retain the existing fixed-dose context and behavior.

GNU Fortran 13.3.0 Debug; Open MPI 4.1.6 for the existing MPI build:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d --parallel 8
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_application_' --output-on-failure
cmake --build build/port-spray-debug --target pelef_spray_3d --parallel 8
python3 tests/spray/test_spray_application.py --exe build/port-spray-debug/pelef_spray_3d \
  SprayApplication.test_piecewise_injection_profile SprayApplication.test_independent_port_rate_profiles
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile tests/spray/test_mpi_spray_application.py
git diff --check
```

MPI: 2/2 application entries pass, 20.64 s, no failures/skips. Added profile
rates [0,0,0.02,0] kg/s at [0,7.5e-9,1.5e-8,2.25e-8] s produce no first-event
dose followed by two 7.5e-11 kg doses. Particle history is 0/4/8/8, final watermark
8 and injected mass 1.5e-10 kg. Injection-corrected mass/energy errors remain below
2e-11. Saved-event continuation on changed ranks (2↔4) matches continuous fields
and particles at relative 2e-11/absolute 1e-25; changed profile input is rejected.
Existing MPI application/schedule/restart cases also pass.

Serial: both existing shared/independent-port profile tests pass (2.615 s), no
failures/skips. Provenance, contract, Python syntax and whitespace checks pass.
No selected-mechanism profile rerun, independent MPI port schedule, full matrix,
physical qualification or production scaling is claimed.
