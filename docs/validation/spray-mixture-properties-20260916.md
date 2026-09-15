# Mixed-liquid component property inputs — 2026-09-16

Tested parent: `65380d5` plus this change. Source/test hashes and build/test/preflight
logs are retained in `build/mixture-properties-checks/`.

## Implementation

The CLI now constructs each mixture component with its own density and heat-capacity
polynomials, bounded Antoine saturation pressure or Clausius surface cap. Existing
component thermodynamics and transport kernels evaluate these models. The effective
injection material uses mixture density, heat capacity and energy evaluated at each
port's temperature, clearing component polynomial and critical-temperature metadata
so the effective properties are not evaluated a second time as one component.

Input validation rejects inactive coefficients, unused component entries, invalid
fits and incompatible vapor-state models. Active models, coefficients and bounds
enter the ordered component restart context; default constant-property contexts stay
unchanged. This adds no physical mixture boiling or wall-film composition model.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d --parallel 8
python3 tests/spray/test_spray_application.py --exe build/port-spray-debug/pelef_spray_3d SprayApplication.test_mixture_property_models_restart SprayApplication.test_mixture_cli_injection_and_restart
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CLI methods passed in 3.733 s, no failures or skips. New cases combine both
component density/heat-capacity polynomials with Antoine or surface-cap settings,
including two injection ports at different temperatures and scheduled injection.
Checks cover analytic mixture diameter at evolved composition/temperature, strict
mass/momentum/energy closure, exact restarted CSV/checkpoint/final-history agreement,
changed-coefficient restart rejection and invalid input rejection. The existing
constant-property initial/shared/independent injection cases also pass. The surface
cap case verifies its CLI/restart connection, not high-temperature cap activation.
Provenance, Python syntax, project contract and whitespace checks passed. No full
suite, configuration matrix or physical validation was run.
