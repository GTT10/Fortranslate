# MPI EB film heat and species exchange — 2026-09-18

Implementation `b3b8528`, base `d08f873`.

The existing distributed film interval now accepts paired EB geometry and
original nodal phi. It migrates films by fluid-side ownership, samples EB gas
state and deposits conserved heat/species using the established EB composite
field API. Adaptive common film substeps and whole-interval rollback are shared
with the regular-wall path. EB mode requires nonperiodic outer faces and a
spatially uniform wall temperature; material conductivity models are retained.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja, existing fixed tests-enabled and selected H2O2 tests-disabled
configurations. MPI test programs use runtime checks and invalid/zero/overflow
floating-point traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_film_step test_spray_mpi_amr_film_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(eb|amr)_film_step_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **4/4 passed**, total **2.26 s**. New EB tests took 0.56/0.58 s
on two/four ranks; regular regressions took 0.56/0.53 s. No test failures or
skips. Selected library compilation, provenance, Python syntax, contract and
whitespace passed. Selected execution was not run.

The first test build failed because the new inventory assertion passed the MPI
communicator in the wrong argument position. Correcting the call to the existing
integration signature fixed compilation; the initial build log is retained.
No production tolerances or golden values changed.

The new fixture has three refinement levels, a 16-cube original nodal field and
a plane at y=0.000375 m. Two synthetic mixed films lie on nested and root-only
surface regions and start on one rank. Both heat and lose mass, with positive
wall heat from a 400 K wall. Composite gas plus film mass, momentum, energy and
species close against wall exchange within 2e-12 times max(abs(initial),1e-12).
Covered gas remains zero and covered temperature retains its sentinel 42.

One allowed substep over a longer interval rejects with unchanged gas, films,
components and ledger. Rank-disagreed EB presence and a nonuniform EB wall
temperature also reject collectively without changes. Empty populations return
zero exchange after geometry/gas validation.

Artifacts: `build/spray-mpi-eb-film-step-checks/{build-initial.log,build.log,build-final.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
Distributed impacts, reactive composition and application persistence remain
separate work. No full matrix, scaling, external parity, installed audit or
physical validation was run.
