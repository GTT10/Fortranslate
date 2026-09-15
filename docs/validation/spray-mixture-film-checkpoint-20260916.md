# Mixed-film component checkpoint — 2026-09-16

Tested parent: `bfedd5f` plus this change. Source/test hashes and build/test/preflight
logs are in `build/mixture-film-checkpoint-checks/`.

## Implementation

Optional `film_component_mass` stores whole-patch component kg, with the same ordered
component table as particle `component_mass`. Film+mixture mode requires both matrices.
The new `_FMIX_V1` suffix distinguishes a component row after each film record.
Shape, finiteness, nonnegativity and summed film mass are validated before writing or
committing a read. Invalid/truncated candidates preserve all external state. Pure-film
and particle-only mixed layouts remain unchanged. The ordered properties and film
settings remain part of caller-managed context. Grid/CLI wiring is still open.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_integration_(wall_coupling|mixture_grid)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two integration entries passed in 0.77 s total, no failures/skips. The wall case adds
exact particle/film composition roundtrip, rejection when the film component output
is omitted, and negative-mass corruption rejection with gas/temperature/component
rollback. Existing pure wall coupling and particle-only mixture checkpoint checks
also pass. The new film checkpoint fixture assigns distinct compositions to existing
film and particle records; it does not prove a mixed-wall CFD trajectory. An initial
regex selected only mixture_grid; the corrected final command above ran both entries.
Provenance, contract and whitespace checks passed; no Python source changed. No CLI,
full configuration matrix or physical validation was run.
