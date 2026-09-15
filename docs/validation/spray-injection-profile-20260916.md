# Tabulated injection rate — 2026-09-16

## Implementation and source

Tested parent: `9bdff9f8b5858db086b2bced963b5ddd9eed30ce` plus this commit's source/tests.
Exact tested hashes and build/test logs are in `build/injection-profile-checks/`.

`references/spray_injection_profile.json` pins the PelePhysics piecewise-linear
ROI interpolation source at `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, with
SHA-256 and license. The application now accepts up to 64 time/rate knots for
its global/shared injection schedule. It integrates the piecewise-linear law
across each batch window, including partial segments, and omits zero-dose
batches while advancing their event counters.

Shared-port fractions divide each integrated dose. Restart reconstructs
cumulative injected mass and the number of nonzero batches, so a pause does
not violate canonical particle/event counts. The complete table enters
checkpoint identity. The no-table constant-rate path retains its context.

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

The CLI CTest entry passed in 40.21 s, including 30 Python tests; no failed or
skipped tests. Build, mechanism provenance, Python syntax and project contract
passed. Compiler warnings include real comparisons.

The new piecewise-linear case contains a triangular pulse, two zero-dose
windows and a rising final pulse. Integrated liquid mass is 4e-11 kg, with
six nonzero batches: 12 parcels for one outlet or 24 for two shared outlets.
Tests check conservation, the analytical mass, omitted zero-dose parcels,
byte-identical gas/parcel CSVs and identical final history after restart from
the pause. Modified rates reject restart; duplicate knots, negative rates and
simultaneous constant-rate input are rejected. Existing constant/independent
schedule tests remain covered.

## Remaining boundary

Each integrated dose is emitted at the beginning of its schedule window;
continuous-time injection trajectories are not claimed. Independent-port
schedules currently retain their constant rates and cannot use this global
table. Per-port rate histories, time-dependent velocity/temperature histories,
physical child creation, mixture droplets and particle MPI/AMR/EB remain open.
No full-suite or external physical validation was run.
