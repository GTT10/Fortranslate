# Mixed-liquid trace and terminal inventories — 2026-09-16

Tested parent: `2186ba4` plus this change. Exact modified source/test hashes and
logs are retained in `build/mixture-terminal-checks/`.

## Implementation

The frozen-gas update transfers an actively evaporating trace component's residual
mass to gas when its current or predicted residual mass is at most 1e-12 of current
droplet mass. Suppressed components remain in liquid. The grid step no longer
uses relative depletion of an already trace component as a timestep restriction.
If the predicted mixture diameter falls below the configured minimum during
evaporation, all component inventories transfer through the existing terminal
energy convention. Evaporation-disabled droplets are preserved.

Pinned Drag.H at `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b` uses a 1e-12
component threshold and a terminal mass cutoff. Its trace path zeroes fractions
and renormalizes. This implementation instead returns the residual species mass
and energy, retaining conservation explicitly. The criterion here uses current
mass for the residual threshold and the configured diameter for terminal removal;
it is not an assertion of exact upstream trajectory parity.

The public evolution tag and CLI restart context change to
`mixed-liquid-components-v2-terminal`. Old mixture checkpoints reject. The file
layout and pure-liquid restart context remain unchanged.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(mixture_terminal|mixture_rates|mixture_exchange)$' --output-on-failure
python3 tests/spray/test_spray_application.py --exe build/port-spray-debug/pelef_spray_3d SprayApplication.test_mixture_terminal_restart SprayApplication.test_mixture_cli_injection_and_restart
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Three CTest entries passed in 0.03 s and the two selected CLI tests passed in
2.549 s, without failures or skips. Provenance, Python syntax, contract and
whitespace checks passed. The full CLI suite was not rerun.

The new cases check trace mass routing and energy closure, ambient-vapor suppression,
complete terminal mass/energy return, evaporation-disabled preservation, and zero
component/diameter CSV fields after full evaporation. Final terminal CSVs and
checkpoints match after restart. A replaced old evolution tag rejects. The existing
mixture injection/restart test still covers initial, shared and independent ports.
No long-time physical validation or full configuration matrix was run.
