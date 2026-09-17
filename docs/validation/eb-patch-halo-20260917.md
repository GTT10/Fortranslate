# EB patch extraction and coarse-donor halos — 2026-09-17

Implementation `2f3dfd2bcad265918e6e4027ab928ef4d9526c58`, based on
`3be2bd4` on local `main`.

## Implemented interface

`extract_eb_geometry_3d` copies an inclusive cell box from valid geometry,
including its surrounding face metrics. Cell spacing and normalized centroid
offsets are unchanged; patch physical bounds are derived from the source box.
Wall centroids remain in global physical coordinates. Invalid boxes or invalid
candidate geometry preserve the existing destination.

`fill_reactive_eb_patch_halo_3d` takes a coarse geometry/state, a fine patch,
its lower coarse-cell index and an integer refinement ratio. Fine dimensions
must be divisible by the ratio and the patch must match the corresponding
coarse Cartesian bounds and spacing. Interior face halos use piecewise constant
parent donors; physical domain face halos extrapolate the nearest fine state.
Only open face slots are populated. Interior, edge, corner and closed-face slots
are zero. Each used donor passes EOS conversion, which provides its temperature.
Nonphysical or solid parent donors and invalid shape/alignment reject without
modifying either halo array.

The API checks grid alignment and donor validity, not complete coarse/fine EB
metric equivalence. Hierarchy integration must ensure consistent restricted
geometry and synchronize fluxes. It currently provides outflow physical faces;
periodic and other physical boundaries, same-level neighbors, temporal donor
interpolation, transport halos and MPI are not implemented by this API.
Allocation failure retains the existing geometry-builder allocation behavior.

## Completed checks

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; existing matching
`build/selected-tree-mpi-debug`, MPI on, spray off. Fixed core and selected EB
runtime compiled; execution uses fixed elementary chemistry. No MPI or application
matrix was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_patch_halo_3d test_eb_geometry_3d test_eb_reactive_transfer_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_patch_halo_3d|unit_eb_geometry_3d|unit_eb_reactive_transfer_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 PASS, 0 failures, 0 skips, 0.21 s**; each test reported 0.06 s.
The new test uses a 12-cubed fine cut-plane geometry, its 6-cubed coarse
restriction and two 6-cubed fine patches: internal and touching domain boundaries.
It checks exact cell/face extraction, every open face donor and temperature,
closed/unused slots, a nontrivial Euler update with the generated halos and
accepted fluxes, and rollback for invalid crop, wrong placement and invalid
parent state. Existing geometry and conserved-state transfer tests also pass.

Provenance and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-patch-halo-checks/{build,ctest}.log`.
The test does not yet evolve a synchronized coarse/fine hierarchy and does not
establish convergence, external PeleC parity or physical validation.
