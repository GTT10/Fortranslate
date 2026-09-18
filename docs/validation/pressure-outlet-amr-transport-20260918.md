# Two-level AMR pressure outlet transport — 2026-09-18

Implementation `31381a2`, base `abc7e41`.

`advance_amr_reactive_transport_3d` now accepts optional `outlet_pressure(6)`
and forwards it to both Euler stages. Coarse physical fluxes and coarse
start/end interpolation halos use the shared pressure ghosts. Fine physical
ghosts apply the target only on touched root faces. The complete two-level
reactive wrapper forwards pressure through both transport halves and hydro,
including its transport-disabled Strang dispatch. Flux accounting, reflux,
average-down, positivity limiting and transactional publication remain intact.

AMR application input and checkpoint identity, patch trees and MPI pressure
connections are still pending. AMR namelist pressure remains explicitly rejected
until its application path is connected.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, Debug/Ninja, existing fixed tests-enabled
and selected H2O2 tests-disabled configurations; no MPI runtime test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_reactive_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_reactive_transport_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The physical test matrix now runs with zero pressure targets and a 102325 Pa
x-high target. It covers an interior patch, face contact, ratio-three corner
and full-domain refinement, with transport-only and inert/reacting full splits.
Existing composite component/Euler and element checks retain their tolerances
(2e-11 component/Euler, 2e-10 elements) and normalization. Invalid wall/inlet
states still must roll back with pressure enabled. Full-domain viscous pressure
ghosts must exactly equal regular-grid pressure ghosts; this is a ghost
comparison, not exact equality of the full AMR transport update.

Artifacts: `build/pressure-outlet-amr-transport-checks/` contains `provenance.log`,
`build.log`, `build-final.log`, `build-selected.log`, and `ctest.log`.

Final CTest **1/1 passed**, 0 skips, 145.60 s total. Fixed and selected builds,
mechanism provenance, Python syntax, project contract and diff checks passed.
The initial fixed build failed at link time because the test lacked the shared
pressure helper import; adding the import fixed it. Tests were launched only
after the final build succeeded. No tolerance changes. This is focused
conservation, rollback and ghost-equivalence evidence, not physical validation.
