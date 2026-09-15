# Mixed-wall CLI integration — 2026-09-16

Tested parent: `d8bd767` plus this change. Source/test hashes and build/test/preflight
logs are in `build/mixture-wall-cli-checks/`; `example/` retains the example outputs.

## Implementation

The CLI accepts mixture wall impact and carries film component matrices through
reactive stepping, total inventories, CSV output and checkpoint read/write. Film
volume checks on restart use current composition and temperature. Film CSV component
columns are kg per whole patch; parcel columns remain kg per physical drop.
Existing wall viscosity/surface tension and conductivity settings apply. Pure-port
property tables remain excluded from mixture mode. The new synthetic mixture_wall
case illustrates deposition and evaporation, not calibrated liquid properties.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d --parallel 8
python3 tests/spray/test_spray_application.py --exe build/port-spray-debug/pelef_spray_3d SprayApplication.test_mixed_wall_cli_restart SprayApplication.test_wall_film_restart SprayApplication.test_mixture_cli_injection_and_restart
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Three CLI methods passed in 9.962 s, no failures/skips. The new method exercises
heated/evaporating deposition, frozen splash and splash with KH children enabled.
It checks strict conservation, film component sum and density/volume, changing film
composition/temperature and exact restarted checkpoint, gas/particle/film CSV and
final-history agreement. Existing pure film and shared/independent mixture injection
checks pass. The standalone `cases/spray_3d/mixture_wall.nml` run completed using
`build/mixture-wall-cli-checks/example/` as working directory, with its log and outputs
retained there. Provenance, Python syntax, project contract and whitespace checks
passed. No full suite, configuration matrix or physical qualification was run.
