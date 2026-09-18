# Two-level AMR pressure outlet hydro — 2026-09-18

Implementation `7ea860b`, base `aa59c20`.

Optional `outlet_pressure(6)` now passes through the serial two-level hydro
and Strang chemistry/hydro APIs. Coarse RK2 uses the regular pressure outlet;
coarse start/end interpolation halos also receive pressure ghosts. Fine RK2
rebuilds pressure ghosts only on physical faces touched by the patch, using
its own adjacent fine state. Interior coarse/fine faces retain interpolation.
Existing subcycle flux accumulation, reflux, average-down and composite
surface accounting remain in place. Updates still publish transactionally.

The shared fine physical ghost builder exposes the same optional argument for
later transport integration. Full transport splitting, AMR application input,
checkpoint signatures, patch trees and MPI pressure dispatch remain pending.
The AMR input reader continues to reject pressure targets until those paths
are connected. This is not full characteristic inflow/outflow or GC-NSCBC.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled
and selected H2O2 tests-disabled configurations. No MPI runtime test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_reactive_chemistry_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_reactive_chemistry_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The added cases use PCM and characteristic PLM with Rusanov, eight patch
placements (six face contacts, full refinement and a ratio-three corner),
and pressure targets on x-high/z-low. They check composite inventory closure
at the existing 2e-11 relative tolerance with floor 1e-20 and common momentum
scale. Full refinement is compared exactly against uniform fine-grid substeps,
including the surface ledger. A nonphysical extrapolated pressure must reject
with unchanged coarse/fine states and temperatures, zero ledger and zero reflux.
The pressure-enabled chemistry split is compiled but not independently run.

Artifacts: `build/pressure-outlet-amr-hydro-checks/`, including `provenance.log`,
`build.log`, `build-final.log`, `build-selected.log`, and `ctest.log`.

Final CTest **1/1 passed**, 0 skips, 69.96 s total. Fixed and selected builds,
provenance, Python syntax, project contract and diff checks passed. The initial
fixed build failed because the added fixture referred to a helper-local config;
using the explicit 102325 Pa test target fixed it. No tolerance changes or
stale-executable runs. Selected build retains existing compiler warnings.
This is focused conservation and uniform-grid equivalence evidence, not
physical validation.
