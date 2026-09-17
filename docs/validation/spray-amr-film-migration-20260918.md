# MPI AMR stationary wall-film migration — 2026-09-18

Implementation: `15c2c8d055dc84f558a3b413543713294e220e59`, based on `2611ee1`.

`migrate_mpi_amr_wall_films` assigns each stationary film patch to the finest
containing cell's owner. Physical bounds are closed for films: exact upper-wall
coordinates choose the final cell without changing stored positions. Exterior
positions reject. Component inventories, geometry, material and thermal state
travel together; source-rank/local-order receive ordering needs no invented film
IDs. There is no patch merging or lateral film transport.

The API collectively verifies topology/owners, bounds, component count and
capacity. Counts precede payload exchange; allocation failures are collective.
A private receive buffer is validated before either caller array is replaced.
Maximum local send/receive count is one million, component count is 1–16. This
is a library ownership step, not integrated AMR wall-film physics or restart.

## Completed checks

GNU Fortran 13.3.0 Debug/Ninja, Open MPI 4.1.6, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`. Reused unchanged fixed tests-enabled and
selected H2/O2 tests-disabled build configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_film --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_film_(2|4)$' --output-on-failure --output-log build/spray-amr-film-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips: two ranks **0.46 s**, four ranks **0.48 s**,
total **0.99 s**. A three-level tree routes lower-wall and two upper-wall
populations, preserves every film/component value exactly, checks deterministic
receive order, changes node owners, coarsens and restores refinement. Component
inventories agree at 2e-12 relative tolerance. Capacity failure, rank-disagreed
owners and exterior-film input reject without committing arrays. Empty ranks
and globally empty film populations succeed.

Fixed and selected builds, provenance, Python syntax, project contract and
whitespace checks passed. No build/test attempt failed. Existing warnings
remain. Selected evidence is compilation only. No Python source was changed.

Logs: `build/spray-amr-film-checks/{build,build-selected,ctest}.log`.
Test executable: `build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_film`.

Remaining integration includes impact-created films, AMR gas sampling and
conservative film evaporation feedback, whole gas/particle/film remap and
application/checkpoint controls. This check does not establish physical impact
validation, EB interaction, convergence or production scaling.
