# Regular 3D characteristic pressure outlet hydro — 2026-09-18

Implementation `b8d2168`, base `6cd6203`.

`advance_reactive_boundary_hydro_3d` accepts optional `outlet_pressure(6)`.
Zero entries keep ordinary ghost filling. Positive finite entries (at most
1e15 Pa) require an existing `outflow` face. Each RK2 stage recomputes both
normal ghost layers from its adjacent physical cell, using distance ratios 2
and 4 in `characteristic_pressure_ghost_3d`. Only face strips needed by the
normal PCM/characteristic-PLM flux stencil are replaced; existing transverse
corner fills remain. Flux-derived surface inventories include the pressure
boundary. Invalid extrapolated states reject before fluid publication.

This connects the pinned nozzle equation kernel, not GC-NSCBC or an acoustic
nonreflection claim. Reaction/transport wrapper, application input, AMR and MPI
pressure-outlet connections are still pending. Reference provenance remains in
`characteristic-pressure-ghost-20260918.md`.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Test serial; selected library built.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_boundaries_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_boundaries_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, 0 skips, 5.84 s total. Builds and static checks
passed. Initial test failed because the added fixture specified unsupported
`plm`; changing it to the existing `characteristic_plm` selection resolved the
failure. No solver tolerance or reconstruction rule changed.

The new cases exercise every face with PCM and characteristic PLM using HLLC,
raising target pressure from 101325 to 102325 Pa. Nonzero boundary exchange is
required and total inventory matches the surface ledger at relative 3e-12 with
unit scale floor. An invalid extrapolated pressure preserves fluid and
 temperature. Existing kernel, reservoir inlet and three-solver boundary cases
remain. This is focused conservation/dispatch evidence, not physical validation.

Artifacts: `build/pressure-outlet-hydro-checks/` contains `provenance.log`,
`build.log`, `build-final.log`, `ctest.log`, `ctest-final.log`, and `build-selected.log`.
