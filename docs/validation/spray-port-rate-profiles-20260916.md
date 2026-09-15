# Independent port rate histories — 2026-09-16

## Implementation and source

Tested parent: `13ac0926abc7b291b32681631067535621dca359` plus this commit's source/tests.
Exact tested hashes and build/test logs are in `build/port-rate-profile-checks/`.

The piecewise-linear ROI source remains pinned in
`references/spray_injection_profile.json`. Independent ports now accept their
own time/rate tables through `port_profile_points`, `port_profile_times` and
`port_profile_rates`. Zero-count ports retain constant rates, so profiles and
constant-rate jets may coexist. Each table must cover its own schedule and
passes finite, bound, strictly increasing time and unused-entry checks.

The existing exact trapezoidal window integrator now accepts a port selection.
Runtime and restart replay use the same merged event order and dose calculation;
zero-dose events advance without creating parcels. Restart counts only
nonzero batches and verifies cumulative injected mass. Tables enter context,
while no-table contexts remain unchanged.

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

The CLI entry passed in 41.08 s, including 31 Python tests; no failures or
skips. Build, mechanism provenance, syntax and contract checks passed. Compiler
warnings include real comparisons.

New cases combine a paused profile with either a constant-rate port or a
second profile, with different injection intervals. Port masses independently
match 4e-11 and 8e-11 kg, and 20 parcels represent the ten nonzero batches.
Checks cover system conservation, byte-identical gas/parcel CSVs and identical
final history after restart from the pause. Modified rates reject restart;
constant/profile conflicts, unused port tables and duplicate times fail.
Existing shared profiles and constant schedules remain covered.

## Remaining boundary

Integrated doses are still emitted at each window's beginning. Tables prescribe
mass rate, not time-varying injection speed, temperature or composition.
Physical child creation/dispersion, mixture droplets, particle MPI/AMR/EB and
physical qualification remain open. No full-suite or external solver run was
performed.
