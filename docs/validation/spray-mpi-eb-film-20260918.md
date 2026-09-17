# MPI EB stationary film ownership — 2026-09-18

Implementation `9bd5ca8`, base `8d5f3a7`.

The existing stationary-film migration API now accepts original nodal EB
geometry. It validates surface location and fluid-facing normal, computes the
finest owner from a fluid-side probe, and transmits unchanged film/component
payloads. The optional argument preserves existing regular-wall behavior.
All ranks agree on EB presence, dimensions and nodal values before exchange.
The nodal lattice may remain finer than the current tree after coarsening.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja. Existing fixed tests-enabled and selected H2O2 tests-disabled
configurations were reused. MPI tests enable runtime checks and floating-point
invalid/zero/overflow traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_film test_spray_mpi_amr_film --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(eb|amr)_film_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **4/4 passed**, total **1.93 s**. New EB tests: two ranks 0.46 s and four
ranks 0.50 s. Existing regular film regressions: 0.46/0.47 s. Selected runtime
library compilation, provenance, Python syntax, contract and whitespace passed.
No failures or skips; selected execution was not run.

The three-level test uses phi=abs(x-2)-1, fluid outside x=1 and x=3, on a 16-cube
nodal lattice. A film exactly at x=1 lies at a nested refinement boundary; the
fluid-side probe assigns it to the nested left patch. Other films exercise the
right patch and unrefined part of the left patch. Changed node owners, root-only
coarsening and subsequent refinement preserve every film field and component
payload exactly. Global component sums retain relative 2e-12 tolerance.

Zero receive capacity, a rank-disagreed owner map, an exterior film, reversed
normal, a rank-disagreed nodal value and rank-disagreed EB option presence all
reject collectively with unchanged payloads. Globally empty films also pass.

Artifacts: `build/spray-mpi-eb-film-checks/{build.log,ctest.log,build-selected.log,provenance.log,contract.log}`.
This is ownership/migration evidence, not distributed EB impact or film-gas
evolution. No application/checkpoint integration, scaling, full matrix,
installed audit, external parity or physical validation was run.
