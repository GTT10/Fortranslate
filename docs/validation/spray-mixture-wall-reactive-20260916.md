# Reactive mixed-wall integration — 2026-09-16

Tested parent: `bb51ace` plus this change. Source/test hashes and build/test/preflight
logs are in `build/mixture-wall-reactive-checks/`.

The reactive driver now accepts film component storage, requires it with mixed film
patches, and carries it through both spray half-steps with gas/particle/film candidates.
Late chemistry or other stage failure commits none of these arrays. Existing pure
property-table and moving impact-wall exclusions remain. CLI wiring remains open.

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_integration_(mixture_wall_coupling|mixture_reactive)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two integration entries passed in 25.96 s total, no failures/skips. New checks cover
reactive mixed wall impact with molecular transport and native implicit chemistry,
committed film composition and total mass/energy with wall budgets. An invalid
chemistry method after the first spray half-step preserves gas, temperatures,
particle/film states and component inventories, with zero returned boundary budget.
Existing mixture reactive checks also pass. Provenance, contract and whitespace
checks passed; no Python source changed. No CLI, new CVODE qualification, full
configuration matrix or physical validation was run.
