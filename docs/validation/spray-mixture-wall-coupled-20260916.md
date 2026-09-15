# Combined mixed-particle wall coupling — 2026-09-16

Tested parent: `3c7a9d6` plus this change. Source/test hashes and build/test/preflight
logs are in `build/mixture-wall-coupled-checks/`.

## Implementation

The mixed particle coupler now accepts wall options, film patches and whole-patch
component masses. It dispatches wall crossing before aerodynamic breakup and skips
breakup on a hit. Retained/splash/deposited inventories preserve composition. New
children and film patches start at the accepted substep end and advance only through
the remaining interval. Particle and film steps share the maximum-substep budget.
Wall impulse, impact energy and film heat enter the boundary ledger. All external
gas/particle/film/component state commits together; capacity or film feedback failure
rolls back. Reactive/CLI integration remains open.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_integration_(mixture_wall_coupling|mixture_grid|mixture_film_grid)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Three integration entries passed in 0.10 s total, no failures/skips in the final
run. New cases cover low-speed deposition with heat/evaporation, high-speed four-child
splash with inherited film composition, total conservation with wall budgets,
remaining-interval film heating and capacity rejection preserving all state. Existing
mixed particle and film grid checks also pass. Initial builds caught a line-length
error and an array declaration preceding its shape argument; both were corrected.
Provenance, contract and whitespace checks passed; no Python source changed. No CLI,
full configuration matrix or physical validation was run.
