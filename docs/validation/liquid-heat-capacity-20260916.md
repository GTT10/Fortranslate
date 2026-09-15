# Liquid polynomial heat capacity — 2026-09-16

Source: parent `2e1ecb1` plus this change. GNU Fortran 13.3 Debug,
`build/port-spray-debug`, fixed H2/O2 mechanism, MPI disabled.

## Implemented scope

Global and per-port cubic heat-capacity fits are validated over the full
specified temperature interval, including interior extrema and the reference
temperature. Energy uses the analytic integral from the existing formation
reference. The frozen-property update and coupling substep estimate evaluate
cp at the initial droplet temperature. This retains the existing approximation
in time; no exact nonlinear thermal integration or new saturation relation is
claimed. Injection and update candidates reject out-of-domain temperatures.
All cp settings participate in restart identity. Constant cp remains the default.

## Focused evidence

Seven CTest entries passed with no failures/skips. Six entries (new heat
capacity, evaporation, coupled integration, rollback, material integration,
application CLI) completed in 20.80 s; the existing constant-cp heat test was
then checked separately in 0.02 s. CLI has 19 tests, with 16.02 s CTest wall time.
The core fixture verifies cp evaluation, the exact energy integral, its numerical
derivative, the frozen-cp heat update and equal/opposite integrated energy
transfer. A quadratic with a negative interior minimum is rejected despite
positive endpoint values. Out-of-range injection is also rejected.
The CLI exercises global and per-port cp fits with active heat and evaporation,
conservation, exact gas/parcel CSV restart and final history equality,
changed-coefficient rejection, invalid fits and unused port inputs.

Build, mechanism provenance, Python syntax, project contract and whitespace
checks passed. No failed runs or broad matrix. These are implementation checks,
not calibrated liquid properties or physical spray qualification.

Artifacts: `build/port-spray-debug/liquid-cp-build.log`, recorded result summary
`build/port-spray-debug/liquid-cp-checks.json`, and the latest constant-heat test
log `build/port-spray-debug/Testing/Temporary/LastTest.log`. CLI fixtures are
temporary. Mechanism source SHA-256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(heat_capacity|heating|evaporation|application_cli|integration_(materials|coupled|rollback))$' --output-on-failure
ctest --test-dir build/port-spray-debug -R '^spray_heat$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The first regex's `heating` alternative matches no test; the actual `heat`
case was therefore run by the second command. Counts above use executed tests.
