# EB PLM twelve-solve CTU flux pipeline — 2026-09-25

Code `95edc3d`, parent `8b8504c`. New `compute_reactive_eb_plm_ctu_fluxes_3d`
connects the three normal PLM Riemann solves, six single-transverse corrected
solves and three final double-transverse corrected solves. Single corrections
use dt/(3 dx); final corrections use dt/(2 dx). Each final transverse direction
uses the corner flux corrected by the other transverse direction. Intermediate
Godunov pressure and velocity feed the pressure-work correction. The staging
was checked against pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, `Source/Godunov.cpp`, read with
`git show` in `/mnt/sn850/Projects/gtt2/external/PeleC`.

The API consumes a working geometry with a two-cell margin around the returned
region, plus the three state/temperature/volume-fraction halo layers required
by the normal predictor. Returned face indices describe working cells 3:N-2,
shifted down by two. Covered cells are skipped, closed normal faces remain zero,
and the single/double aperture gates are retained. Optional primitive and
internal-energy sources are applied at half-step only in the final stage.
All stages use private storage; final states, energies, fluxes and face data
publish together. Temporary storage currently retains full normal/corner fields.

These are intensive face-center fluxes. Artificial viscosity, centroid
interpolation, conservative time update/redistribution and application/MPI
halo preparation/dispatch remain. Flattening is still caller supplied. This
increment does not claim a complete selectable EB Godunov solver or external
full-solver parity.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, reused Ninja Debug
`build/spray-amr-reactive-debug`. **1/1 CTest passed**, 1.86 s (1.90 s total),
no failures/skips. The EB regression now checks a regular linear density field
advected by a three-component velocity: all directional mass fluxes at an
interior face agree with the unsplit half-time analytic value within 2e-12.
It also checks x/y/z cut geometries with uniform flow against analytic Euler
fluxes within 3e-13, cropped output bounds, closed-face zeros, exact zero-source
parity and failed-stage rollback. Existing normal/transverse/EB hydro checks
remain. No existing tolerance or golden was changed.

Provenance, Python syntax, project contract and whitespace checks passed.
Compiler warnings remain. No full configuration matrix, selected mechanism,
application/MPI CTU, compiled external oracle or physical validation was run.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-ctu-flux-checks/{build.log,LastTest.log}`.
