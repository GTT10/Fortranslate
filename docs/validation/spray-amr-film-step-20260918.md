# MPI AMR mixed-film heating and evaporation — 2026-09-18

Implementation: `7c641229d3851c5049fef59a6eef72b7a05c4ffb`, based on `9f1d0fe`.

`advance_mpi_amr_film_coupling` connects existing mixture-film rates/exchange
to owner-local AMR gas sampling, conservative deposition and EOS recovery.
It first privately migrates films to their finest owners. Axis-aligned film
normals must select enabled physical walls, with exact wall-normal positions.
Gas query points clip to half the global finest spacing inside the root domain;
physical film positions are not moved. The existing composite stencil remains
piecewise constant across refinement interfaces.

Global minimum substeps limit component loss and temperature change. Rejected
film/feedback trials halve up to 30 attempts. Gas, temperatures, film/component
arrays and signed wall-energy ledger publish together after the whole interval
and final migration succeed. Incoming state and ledger survive failure.
Tables/options are shared caller contracts; face agreement, interval, component
count and evolution mode are checked collectively. Local sample count is bounded
to 125000, including after initial film migration.

## Completed focused verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`. Reused unchanged fixed tests-enabled and
selected H2/O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_film_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_film_step_(2|4)$' --output-on-failure --output-log build/spray-amr-film-step-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips: two ranks **0.52 s**, four ranks **0.53 s**,
total **1.10 s**. The three-level fixture has mixed water/peroxide films at both
x walls, initially all on the last rank, with 1000 K gas, 330 K films and 400 K
walls. Films heat and lose mass; exact physical wall positions remain unchanged.
Independent composite gas-plus-film inventories close mass, momentum, species
and energy including wall heat within 2e-12 of `max(1e-12,abs(initial))` per
component. A long interval with one allowed substep rejects and preserves gas,
temperature, film geometry/state/composition and sentinel ledger. Globally empty
films return zero exchange successfully.

Fixed and selected builds passed; selected evidence is compilation only.
Provenance, Python syntax, project contract and whitespace checks passed.
No build/test attempts failed. Existing compiler warnings remain; no Python
source changed.

Artifacts: `build/spray-amr-film-step-checks/{build,build-selected,ctest}.log`.
Executable: `build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_film_step`.

This is a film exchange interval. It does not yet connect impact-created films,
reacting gas splitting, combined remap, application controls or film checkpoint
state. EB impact, lateral film transport, convergence, external solver parity,
physical validation and production scaling remain open.
