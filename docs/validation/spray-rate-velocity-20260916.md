# Injection velocity from rate — 2026-09-16

## Source and implementation

Tested parent: `a24f33c324e6205ca76d643c74bb67bb2eb3172b` plus this commit's source/tests.
Exact tested hashes and build/test logs are in `build/rate-velocity-checks/`.

`references/spray_injection_profile.json` now pins SprayJet.cpp's readROI
velocity relation at PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
Optional `velocity_from_rate` computes U=mdot/(rho*A*Cd), using each batch's
average flow, its actual nozzle area and the appropriate liquid density at
injection temperature. Global and per-port discharge coefficients are bounded
and enter checkpoint identity. The existing prescribed-speed mode remains
the default.

Single, shared and independent schedules use the same relation. Shared flow
is divided by port fractions; independent flow uses each port interval.
Peak rates are checked through the injection preflight. Zero radius, invalid
Cd and initial liquid inventory are rejected for the new scheduled-only mode.

## Focused verification

GNU Fortran 13.3.0, existing Debug spray build, MPI disabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_application_cli$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The CLI entry passed in 42.96 s, including 32 Python tests; no failures or
skips. Build, mechanism provenance, syntax and contract checks passed.
Compiler warnings include real comparisons.

The new cases check single, shared two-port and independent two-port rates,
with different port densities, radii and discharge coefficients. Every parcel
speed agrees with independently reconstructed batch rate divided by rho*A*Cd
(relative tolerance 2e-13), and speeds vary with rate. System conservation,
byte-identical gas/parcel CSVs and identical final history after restart are
checked. Changed coefficient rejects restart; point outlet, zero coefficient
and initial inventory are rejected.

## Remaining boundary

One average speed and a beginning-of-window emission represent each batch.
This is not continuous mass-weighted momentum-flux integration or external
trajectory parity. Explicit time-dependent temperature/composition histories,
physical child generation/dispersion, mixture droplets and particle MPI/AMR/EB
remain open. No full-suite or physical validation was run.
