# Mixed-liquid film exchange — 2026-09-16

Tested parent: `df95899` plus this change. Source/test hashes and build/test/preflight
logs are in `build/mixture-film-exchange-checks/`.

## Implementation

The stationary-film exchange API accepts component mass losses, vapor enthalpies,
gas heat and wall heat. It solves the remaining mixture energy for temperature and
updates height using mixture density at fixed footprint. All component masses are
whole-patch kg. Gas species gains equal liquid component losses; gas energy plus
wall energy balances the change in film energy. Partial, complete and minimum-height
terminal removal use explicit inventory accounting. Invalid mass loss or failed
property/temperature evaluation preserves caller state and returns zero budgets.

Gas/wall transfer-rate evaluation and the mixed-film grid/checkpoint path are still
open. This is not physical film validation or a completed mixed wall application.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(mixture_wall_exchange|film_evolution)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CTest entries passed in 0.03 s total, no failures or skips. Mixed checks cover
unequal component losses, gas species gains, combined film/gas/wall energy closure,
density-based geometry, overdraw rollback, full evaporation and minimum-height
removal. Existing pure-film evolution checks also pass. Provenance, project contract
and whitespace checks passed; no Python source changed. No CLI, full configuration
matrix or physical validation was run.
