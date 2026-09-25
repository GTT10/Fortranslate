# Single/double EB transverse state correction — 2026-09-25

Code `6c06316`, parent `2665369`. New `apply_reactive_eb_transverse_3d`
ports the cell-local `pc_transdo` / `pc_transdd` correction and species updates
from pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`Source/Godunov.H:30-134,455-849`, inspected with `git show` in
`/mnt/sn850/Projects/gtt2/external/PeleC`.

Inputs are two predicted normal primitive states and internal-energy densities,
one/two pairs of transverse intensive fluxes, Godunov pressure/velocity,
coefficients, gamma and normal/transverse apertures. The single correction
requires both transverse faces open; the double correction requires all four
open. If the joint gate fails, transverse data are not read. Apertures are
Boolean gates here, not weights. Closed normal states are retained.

The correction updates density, momentum, total energy and species densities,
then recovers velocity and internal-energy density. Pressure uses the separate
Godunov pressure-work formula, with ordered subtraction of the two directions.
Double correction can add a complete optional primitive-source group
(source vector, internal-energy source, source interval). Density/species use
the pre-source corrected density for division, matching the reference order.
Invalid inputs/nonpositive corrected density reject without publishing either
normal state or internal energy. Species are not renormalized in this kernel.

This is a cell-local API; the full EB CTU orchestration, corrected Riemann solves,
centroid interpolation and application time stepping remain to be connected.
No full Godunov implementation or external numerical parity is claimed here.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, reused Ninja Debug
`build/spray-amr-reactive-debug`. **1/1 CTest passed**, 1.04 s (1.08 s total),
no failures/skips. The existing EB hydro regression adds every transverse
open/closed mask for one/two directions crossed with four normal masks (80
combinations). Inactive transverse arrays contain NaNs. Checks cover the joint
gate, normal-face retention, density/momentum/total-energy/species conservation,
analytic ordered pressure work, primitive source addition and density-failure
rollback. Existing EB normal-prediction, flux and hydro checks remain. No
existing tolerances or golden outputs changed.

Provenance, Python syntax, project contract and whitespace checks passed.
Compiler warnings remain. Full configuration matrix, application/MPI CTU,
compiled external kernel oracle and physical validation were not run.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-transverse-checks/{build.log,LastTest.log}`.
