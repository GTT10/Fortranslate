# Exact nodal-face symmetry and recursive MOL restart — 2026-09-27

Implementation `19d167b8` preserves an exact geometric invariant. If both pairs
of face endpoint level-set values are identical along a transverse direction,
the clipped region is extruded along that direction and its centered first
moment is exactly zero. The builder now assigns that exact zero instead of
retaining cancellation error from summing two clipped-triangle moments. The
predicate uses exact endpoint equality, with no geometric tolerance. General
oblique faces retain their original calculation. Areas and cell volumes are
unchanged.

## Cause established in the failing fixture

Temporary diagnostics in `build/mol-negative-origin/case/run1..4.log` (the
first file is named `run.log`) isolated the earlier checkpoint rejection:

- At fine cell `(2,5,3)`, old fuel density and old fuel rate were zero.
- New-stage diffusion and particle rates were zero. The MOL hydro rate was
  `-1.9662077639170553e-35 kg/m³/s`.
- All six raw face fuel fluxes were zero. Centroid interpolation introduced
  outgoing z-face fuel flux `9.8310388195852779e-39 kg/m²/s`.
- The trapezoidal update over `5e-8 s` produced
  `-4.9155194097926379e-43 kg/m³`. Coarse averaging carried it into the root
  patch; particle source correction had not yet run.

The exactly x-invariant-in-y/z planar fixture should not acquire transverse
centroid displacement from arithmetic cancellation. Enforcing the exact
symmetry removes this artificial flux and resolves the same fixture's
checkpoint failure. Diagnostic instrumentation was removed before the final
build. No cell species clipping, checkpoint tolerance change or flux acceptance
relaxation was made. This repair does not prove positivity of arbitrary
cut-cell MOL updates.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing selected heptane Simple-transport
Debug MPI/spray build, tests/SUNDIALS disabled:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_core pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2 > build/mol-negative-origin/repair-build.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/gfortran -O0 -g -fcheck=all -ffree-line-length-none -I build/heptane-simple-mpi-debug/modules tests/unit/test_eb_geometry_3d.F90 build/heptane-simple-mpi-debug/libpelef_core.a -o build/mol-negative-origin/geometry-check-standard
build/mol-negative-origin/geometry-check-standard > build/mol-negative-origin/geometry-final.log 2>&1
python3 build/mol-negative-origin/run_repaired.py > build/mol-negative-origin/repaired.log 2>&1
python3 -m py_compile build/mol-negative-origin/run_repaired.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All final commands passed. The existing geometry executable includes six new
checks: valid nodal planes and exactly zero invariant-direction face centroids
for each of three axes. Existing oblique, complementary-volume and coarsening
checks also pass. Its first manual run enabled invalid-operation trapping and
stopped at the existing intentionally NaN-valued invalid-geometry test. The
standard run without that extra trap passes; the failed log remains.

Six nested recursive application launches pass: continuous MOL, first-step
checkpoint, one-to-two-rank restart, MOL-to-CTU rejection, CTU checkpoint, and
CTU-to-MOL rejection. The original partial-refinement 4³ root/ratio-two child
case retains four evaporating parcels, two `1e-7 s` root steps and wall heat.
Gas CSV bytes and sorted particle rows match across restart. Maximum absolute
combined component residual is `2.710505431213761e-20` (component units differ).
Commands and budgets are in
`build/mpi-mol-recursive-particle-symmetric/result.json`; binary SHA256 is
`160bc2ea44c7a045d1a1004b314e57cc3c86fc783b4bb8a6fed4449f53b4175f`.

This supersedes the specific failure in the
[recursive-particle diagnostic](mpi-mol-recursive-particle-20260927.md).
Reacting/film MOL, full-duration reference cases, arbitrary-geometry positivity,
the independent Jet failure and full upstream temporal parity remain open.
