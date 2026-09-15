# Mixed-liquid CLI and injection — 2026-09-16

Tested parent: `85542b4` plus this change. Exact modified source/test hashes,
build/test logs and the example's outputs are retained in `build/mixture-cli-checks/`.

## Implementation

The namelist accepts a shared component table, initial mass fractions, constant
component density/cp and reference Clausius/energy properties. Components and
fractions are validated before outputs; unused entries, duplicate gas species
and invalid properties reject. Initial and scheduled injections create the
component-mass matrix using the requested mixture density and energy at each
port temperature. Per-port identity continues to identify the injection port.

Shared and independent port schedules append component masses and species-resolved
injection ledgers. Rate-derived speed uses the mixture injection density. The
application dispatches mixed parcels through the reacting-flow mixture path,
computes mixture-aware total inventories/diameters, emits component mass CSV
columns, and writes/reads the composition checkpoint extension. Ordered component
properties and fractions participate in restart compatibility. Pure-liquid
checkpoint arguments retain their absent/default behavior.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d --parallel 8
python3 tests/spray/test_spray_application.py --exe build/port-spray-debug/pelef_spray_3d SprayApplication.test_mixture_cli_injection_and_restart
ctest --test-dir build/port-spray-debug -R '^spray_application_cli$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The full CLI entry passed all 35 Python tests in 50.05 s (50.07 s CTest total),
without failures or skips. Provenance, Python syntax, contract and whitespace
checks passed. Development compilation first needed a long namelist continuation
split; the initial port fixture also omitted required nozzle axes, which was
corrected without relaxing input validation.

The new test covers initial mixed droplets, shared two-port schedules with distinct
injection temperatures, and independent two-port schedules with discharge-derived
speed. It checks component-mass closure, diameter from current composition,
preferential evaporation, injection-dose accounting and system conservation.
Continuous/restarted gas and parcel CSVs and final checkpoints are byte-identical;
final histories match. Changed composition rejects restart; malformed fractions,
densities, duplicate species and unused entries reject.

The checked-in `cases/spray_3d/mixture.nml` was also executed from the artifact
directory. It completed 12 steps to 1 microsecond. Maximum absolute mass and energy
errors were 1.42e-23 kg and 3.20e-18 J; the largest momentum error was
3.49e-25 kg m/s. Its H2O/H2O2 properties are synthetic implementation fixtures,
not a calibrated liquid-mixture dataset.

## Remaining work

The composition table is shared by all ports. Per-port mixture compositions,
component property-fit/Antoine/boiling-cap CLI inputs, mixture wall impact and KH
shed children remain open. Representative ETAB/KHRT and computational splitting
are available with prescribed mixture viscosity/surface tension; the preceding
API tests supply their focused evidence. No particle MPI/AMR/EB, full configuration
matrix, new CVODE qualification or physical validation is claimed.
