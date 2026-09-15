# Per-port liquid compositions — 2026-09-16

Tested parent: `3f1e261` plus this change. Source/test hashes and build/test/preflight
logs are retained in `build/port-mixture-checks/`.

## Implementation

`port_mixture_fractions(component,port)` overrides the shared injection composition
for individual ports. Zero columns inherit the shared composition; nonzero columns
must be nonnegative and normalized, with inactive rows/columns zero. The ordered
component property table remains common to all ports.

Initial and scheduled particle inventories, effective injection density/energy,
rate-derived injection speeds and injected component budgets now use the selected
port composition. Existing material IDs identify ports. Active compositions enter
restart identity; default shared-composition contexts remain unchanged. Subsequent
transport evolves each particle's stored component masses independently.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d --parallel 8
python3 tests/spray/test_spray_application.py --exe build/port-spray-debug/pelef_spray_3d SprayApplication.test_port_mixture_compositions_restart SprayApplication.test_mixture_cli_injection_and_restart
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final run: two CLI methods passed in 3.853 s, no skips. New cases cover explicit
20/80 composition, inherited 60/40 composition, a pure component port, initial and
shared/independent scheduled injection, and rate-derived speed. With evaporation
disabled, every output particle retains its expected fraction and analytic mixture
diameter. Mass/momentum/energy checks and exact restarted checkpoint/CSV/final
history comparisons pass. Changed composition and invalid inputs reject. The
existing shared-composition evaporation/injection/restart cases also pass.

Two initial test attempts failed on invalid test inputs: rate-derived speed with
initial droplets, then independent schedules with nonzero global flow. The fixtures
were corrected to scheduled-only injection and zero global flow, preserving those
existing guards. Provenance, syntax, contract and whitespace checks passed. No full
suite, configuration matrix or physical validation was run.
