# Mixed-liquid film rates — 2026-09-16

Tested parent: `c81fdf1` plus this change. Source/test hashes and build/test/preflight
logs are in `build/mixture-film-rates-checks/`.

## Implementation

The mixture skin builder now optionally uses the existing wall half-composition and
half-positive-temperature-difference convention. The default droplet rule is unchanged.
`mixture_wall_film_rates` feeds that skin into the existing stationary-film transfer
correlation and partitions evaporation with component diffusion weights. Its thermal
rate uses component-specific vapor-minus-liquid enthalpies. Liquid conductivity is
supplied, so calibrated mixture conductivity remains a caller responsibility.

`advance_mixture_wall_film_frozen` evaluates rates once and applies their integrated
heat/component transfers with the conservative film exchange API. Component overdraw
or invalid temperature rejects the interval without changing state. The driver must
reduce the step; no component clipping hides overdraw. Grid/CLI/crossing/restart
connections remain open. These are bounded APIs, not physical spray qualification.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(mixture_wall_rates|mixture_rates|film_evolution)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Three CTest entries passed in 0.03 s total, no failures or skips. New checks cover
positive two-component heat/evaporation rates, component-specific thermal rate,
frozen component losses and gas/film/wall energy closure, excessive-step rollback,
and agreement with pure wall-film mass/heat rates when the second component is
absent. Existing mixture droplet rates and pure film evolution also pass. Provenance,
contract and whitespace checks passed; no Python source changed. No CLI, complete
configuration matrix or physical validation was run.
