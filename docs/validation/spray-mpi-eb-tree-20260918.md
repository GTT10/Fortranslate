# MPI EB tree particle ownership — 2026-09-18

Implementation: `334a646e73af2166256be286b2f332876dee25c7`, based on `f7b7805`.

`migrate_eb_tree_spray_parcels` wraps the existing AMR particle migration with
private population/component arrays. It requires a replicated original nodal
field spanning the root domain at the finest tree resolution, validates exact
field agreement across ranks, and checks received point locations with the
tetrahedral level-set query before committing on any rank. Regular topology,
owner-map, domain, payload and capacity checks remain active. Physical outer
ownership is half-open and nonperiodic. The caller binds the field to the gas EB
hierarchy; rounded/coarsened gas volume metrics are not independently checked.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; unchanged existing fixed
tests-enabled and selected H2/O2 tests-disabled configurations reused.
Compiler `/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`, launcher
`/usr/bin/mpiexec`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_tree --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_tree_[24]$' --output-on-failure --output-log build/spray-mpi-eb-tree-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, two ranks **0.48 s**, four ranks **0.49 s**, total **1.02 s**.
Selected library compilation, provenance, syntax, contract and whitespace
checks passed. No build/test failures or skips; existing real-comparison warnings
remain. Selected execution and a broad MPI matrix were not run.

The fixture has a 4-cubed root with one ratio-two child, a replicated 8-cubed
nodal geometry and planar EB x=0.75. It checks finest-child ownership, movement
to an exposed coarse cell, reassigned gas owners, complete ID/velocity/component
values and an empty population. A rank-local geometry perturbation, one solid
particle and zero receive capacity each reject collectively while preserving
the original local IDs, positions and component masses. The solid-point failure
occurs after private communication, exercising whole-operation rollback.

Artifacts: `build/spray-mpi-eb-tree-checks/{build.log,build-selected.log,ctest.log}`.

This is particle ownership groundwork, not distributed EB sampling/deposition,
particle time stepping, gas coupling or rank-neutral EB spray restart. The
replicated nodal geometry and full-field agreement are bounded to one million
finest cells and have no production scaling qualification. No convergence,
external parity or physical validation was performed.
