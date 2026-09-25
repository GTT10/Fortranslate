# Physical boundary halos for EB CTU — 2026-09-25

Source base `f7d1ff9`; the implementation commit immediately precedes this
record. PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`advance_reactive_eb_boundary_ctu_3d` prepares six state-halo layers for explicit
six-face boundary labels. It reuses the regular 3D donor mapper for periodic
wrapping, fixed-wall normal reflection, prescribed inflow and outflow. It forwards
the labels to the existing physical Riemann correction and no-slip viscosity
mask. The outflow entry delegates with six outflow labels, preserving its API.
This is the hyperbolic hydro boundary operation; an `isothermal_wall` label here
reflects normal velocity and does not impose a thermal transport temperature.

The caller supplies compatible target-resolution extended EB geometry, including
periodic metrics where relevant, and prescribed inflow U/T for active inflow
faces. Open target cells requiring a covered state donor reject. All state and
ledger outputs remain private until success. No cell-composition checks, NASA
bounds or tolerances were relaxed.

## Focused evidence

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing
`build/spray-amr-reactive-debug`. `regression_reactive_eb_hydro_3d` passed **1/1**
in 24.17 s (24.22 s total), no runtime failures/skips/interrupted runs.

The existing EB suite now additionally checks five boundary labels in all three
axes: periodic, slip wall, no-slip wall, isothermal wall and specified inflow.
Each checks composite inventory closure at the unchanged `2e-12` scaled criterion;
all fixed-wall normal mass fluxes are exactly zero. Three cut-plane orientations
also pass with tangential physical walls. The six-layer width exceeds the four
physical cells, exercising repeated reflection/wrapping. Missing inflow data and
unpaired periodic labels reject without publishing output. Existing outflow
redistribution, prepared CTU, flattening and physical Riemann tests also pass.
The initial build rejected one new overlong declaration; corrected before testing.
Compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace checks passed. Mechanism SHA unchanged:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-ctu-boundary-halo-checks/{build.log,build-final.log,build-tests.log,LastTest.log}`.

This adds a library boundary-preparation entry. Coupled transport boundary input,
application selection, tree/MPI general physical boundary preparation and moving
walls remain separate work. No full reference or physical wall-flow qualification.
