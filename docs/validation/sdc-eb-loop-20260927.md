# Coupled single-grid EB SDC — 2026-09-27

Implementation `ec33d58b`, base `a5a8fb02`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

`advance_reactive_eb_sdc_3d` couples old/new molecular transport and external
forcing, source-aware EB CTU hydro and forced chemistry from the same old
interval state. Each iteration predicts with averaged nonreacting sources plus
persistent reaction estimates, evaluates new sources, then integrates reaction
with averaged nonreacting rates plus hydro. Hydro and transport source stages
use the selected existing redistribution algorithm before entering the loop.
`none`, `flux_redist` and `state_redist` are supported.

The physical geometry must match the central region of the extended geometry.
Hydro boundary filling and transport physical ghosts share supplied boundary
labels and inflow data. Transport wall temperature is forwarded. Covered data
are excluded from arithmetic and restored at publication. State, temperature,
reaction estimates, integrated fluxes, limiter and all inventories remain
private until the complete interval succeeds.

The exported flux is hydro plus trapezoidal old/new transport. Domain boundary,
EB wall and external-source inventories are separate. Transport boundary
inventory is independently summed from open-area integrated face fluxes;
external work uses trapezoidal old/new momentum-derived velocity.

## Focused evidence

Artifacts: `build/sdc-eb-loop-checks/` (local). GNU Fortran 13.3.0, existing
Ninja Debug tests-enabled build `build/spray-amr-reactive-debug`;
compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_reactive_eb_chemistry_3d
  --parallel 2` (`build.log`). Compiler warnings remain.
- Final CTest **1/1 passed**, 3.15 s (3.20 s total): `env -u LD_LIBRARY_PATH
  -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R
  '^unit_reactive_eb_chemistry_3d$' --output-on-failure` (`ctest.log`).
- Added a plane-cut 2-cubed physical fixture with extended geometry, two SDC
  iterations, dt=1e-9 s, elementary H2/O2, all transport flags, forcing
  `[20000,-30000,40000]` and outflow faces. All three redistribution modes
  complete. The first five conserved-component changes match domain boundary,
  wall and forcing inventories within 2e-10 times max(1, initial absolute
  inventory); active density equals summed species within 2e-12.
- Invalid artificial viscosity after old transport evaluation rejects the
  whole interval, preserving state, temperature, reaction estimates, fluxes
  and sentinel diagnostics exactly. Existing masked chemistry, forced reaction
  and Strang checks also passed.
- An earlier run passed before transport boundary accounting was changed
  from rate integrals to independent face-flux sums; the final run above is
  the reported evidence. Python syntax for `tools/check_project_contract.py`,
  project contract, touched Fortran line lengths and `git diff --check` passed.

This is bounded conservation/rollback evidence, not pinned whole-solver parity
or a temporal-convergence qualification. No application/restart, disabled-
operator coupled case, other coupled boundary case, AMR/MPI EB SDC, selected
build, full matrix, installed audit or physical validation was run. The full
port remains incomplete.
