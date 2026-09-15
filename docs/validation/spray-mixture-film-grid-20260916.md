# Mixed-film grid feedback — 2026-09-16

Tested parent: `177b6d1` plus this change. Source/test hashes and build/test/preflight
logs are in `build/mixture-film-grid-checks/`.

## Implementation

The stationary mixed-film grid API samples gas at each patch, computes heat and
component evaporation rates, limits depletion/temperature substeps and scatters
component and energy increments to gas. Failed gas recovery retries with half the
interval; ultimate failure preserves every external array and returns zero budgets.
Optional per-patch start times support film births during a later particle step.
Terminal patches and their component columns are removed together. The totals API
now includes optional stationary film species and internal energy. Particle/reactive
and CLI integration remains open; this does not establish a complete mixed-wall run.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_integration_(mixture_film_grid|mixture_grid)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final run: two integration entries passed in 0.08 s total, no skips. New checks
cover evolving mixed films, mass/species/energy closure with wall heat, substep-limit
rollback, rejection of terminal release exceeding available gas energy, and removal
of an admissible terminal film with its component column. Existing mixed particle
coupling checks also pass.

During implementation, a missing test import caused a link failure; allocating the
candidate temperature before the recovery call fixed an initial-state rejection.
The large-film terminal fixture correctly failed gas admissibility, so it became an
explicit rollback case and a smaller film separately exercises successful removal.
Tolerances and gas recovery guards were retained. Provenance, project contract and
whitespace checks passed; no Python source changed. No CLI, full configuration matrix
or physical validation was run.
