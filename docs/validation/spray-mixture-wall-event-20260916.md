# Mixed-liquid wall-impact event — 2026-09-16

Tested parent: `bbfd884` plus this change. Source/test hashes and build/test/preflight
logs are in `build/mixture-wall-event-checks/`.

## Implementation

`resolve_mixture_wall_impact` uses the existing six-regime event kernel with mixture
density, heat capacity and energy evaluated at the incident temperature. Retained
parent and splash children inherit per-drop component masses; the deposited film
receives total component masses for its complete patch. Temperature and composition
are unchanged across the instantaneous event. Wall impulse/energy account for
velocity changes and stationary deposition. Failed construction commits no particle,
film or component state. No film entrainment is introduced.

This is an event API with caller-supplied dimensionless parameters, not the completed
mixture wall path. Boiling-temperature evaluation, planar crossing dispatch, evolving
film composition and grid/CLI/checkpoint connection remain open. Existing impact
correlation admissibility limits remain unchanged.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(mixture_wall_event|wall_event)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CTest entries passed, no failures or skips (0.03 s total). New mixed cases cover
all six regimes with a component heat-capacity polynomial, component conservation,
inherited splash composition, total energy and wall momentum closure. A late
film-construction failure preserves parent/child/film component inventories and
particle/film state. Existing pure event cases also pass. Provenance, project contract
and whitespace checks passed; no Python source changed. No CLI, full configuration
matrix or physical impact validation was run.
