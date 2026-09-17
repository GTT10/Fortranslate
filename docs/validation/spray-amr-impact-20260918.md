# MPI AMR planar impact and film creation — 2026-09-18

Implementation: `58eda2cd72cd51834f845665978246d339f38d38`, based on `64fd6e1`.

The AMR particle interval accepts optional film arrays/settings together with
explicit boundary faces and an ID watermark. It calls the existing distributed
planar impact model after parcel/breakup updates and before generic boundary
handling. KH children inherit parent pre-step locations and sampled pressure
for impact queries. Existing thickness queries, impact classification, deposition
and optional splash creation retain their existing model semantics.

Trial films remain private and are accepted only with particle/gas feedback.
Dead parcels are removed between accepted impact substeps. Film state migrates
to finest AMR owners before final publication. Global wall impulse/energy is
accumulated separately from rank-local flow-boundary exchange, then subtracted
from the final system ledger without another rank sum. Whole-interval failure
preserves gas, particle/film components and IDs. The reacting wrapper passes its
private films through both particle halves and now requires `last_id` in film
mode, including when the initial airborne population is empty.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`. Reused unchanged fixed tests-enabled and
selected H2/O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_impact test_spray_mpi_amr_film_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(impact_(2|4)|film_reactive_2)$' --output-on-failure --output-log build/spray-amr-impact-checks/ctest.log
# After correcting the new impact fixture:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_impact_(2|4)$' --output-on-failure --output-log build/spray-amr-impact-checks/ctest-impact-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final impact tests **2/2 passed**, two/four ranks **0.55/0.55 s**, total **1.14 s**.
The existing two-rank film/reacting-wrapper composition regression passed
**20.89 s** in the earlier combined invocation. No skips. Final fixed and
selected builds, provenance, syntax, contract and whitespace checks passed.
Selected evidence is compilation only; existing warnings remain.

The three-level impact fixture disables drag/heating/evaporation to isolate
boundary accounting. One parcel exits the high-x outflow, one deposits a
represented mass of 1e-10 kg on the low-x wall, and one remains airborne.
Independent gas-plus-particle-plus-film mass, momentum, energy and species
changes match the boundary ledger within 3e-12 of inventory/ledger scale.
Assertions check one film, its exact wall position and represented mass, and
surviving parcel identity. Substep exhaustion, rank-disagreed face settings and
insufficient parcel capacity preserve the original film-free state and ledger.
The earlier composition regression exercises film evaporation with gas transport.

Initial failures: the new test declared film and AMR patch arrays with the same
name, causing compilation failure; they now have distinct names. The first
runtime attempt omitted nonzero surface tension and liquid viscosity required
by the impact model, so the impact step rejected. Explicit synthetic values
0.072 N/m and 0.001 Pa s fixed the fixture. No tolerance/model change was used
to turn those checks green.

Artifacts: `build/spray-amr-impact-checks/` contains `build-test.log` (initial
compile failure), `build-final.log`, `build-impact-final.log`, `build-selected.log`,
`ctest.log` (composition pass plus initial impact failures) and
`ctest-impact-final.log`. New executable:
`build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_impact`.

This newly tests deposition, not every rebound/splash regime or simultaneous
KH/splash creation. Film-aware application controls and checkpoint state remain
open. EB impacts, lateral film transport, convergence, experimental validation,
external parity and production scaling are not established by these checks.
