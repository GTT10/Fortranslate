# EB face-centroid flux interpolation — 2026-09-25

Code `0d110d5`, parent `2aa0b5d`. New
`interpolate_reactive_eb_centroid_fluxes_3d` implements the 3D bilinear
face-center-to-centroid interpolation in pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, `Source/Hydro.H`
(`eb_compute_div`, approximately lines 314-432). Source was inspected with
`git show` in `/mnt/sn850/Projects/gtt2/external/PeleC`.

Inputs include one surrounding cell of geometry and face-center intensive
fluxes. Outputs correspond to central cells 2:N-1, shifted down by one. For a
partially open face, the signs of the two centroid offsets select tangential
neighbors and the absolute offsets give interpolation weights. Axis-neighbor
connectivity is obtained from open tangential faces of the positive-side cell,
including the upper output boundary's halo cell. Closed/zero-weight donors
are skipped. Regular faces use their center value and closed faces return zero.
Outputs publish together only on success.

Tracing `pc_adjust_fluxes_eb` through `pc_ext_flx_eb` confirmed that the reference
multiplies by regular face area only before interpolation; it does not multiply
by aperture there. On a Cartesian grid the regular area is constant along each
face direction, so this API interpolates intensive flux directly. The existing
EB divergence then multiplies by target aperture exactly once. The initial
interpretation that aperture-weighted flux should be interpolated was corrected
before building/testing or committing.

The API is not yet connected to the new CTU pipeline or time update. Viscosity
must precede this interpolation when that connection is made. Exact AMReX
connectivity parity for more complicated topology remains unqualified.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, reused Ninja Debug
`build/spray-amr-reactive-debug`. **1/1 CTest passed**, 0.07 s (0.11 s total),
no failures/skips. The EB wall-flux unit test adds regular geometry and two
oppositely oriented oblique planes, all three face directions, constant and
bilinear tangential fields. Interpolated values match evaluation at the actual
centroid within 1e-14 (constant) / 1e-12 (bilinear); both nonzero tangential
offsets are required. Closed donor data are NaNs. A bad active donor preserves
all three output arrays. Existing wall-flux/geometry checks remain. No existing
tolerance or golden changed.

Provenance, Python syntax, project contract and whitespace checks passed.
Compiler warnings remain. CTU-to-centroid integration, disconnected-neighbor
runtime fixtures, application/MPI execution, compiled external oracle and
physical validation were not run.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_reactive_wall_flux_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_reactive_wall_flux_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-centroid-flux-checks/{build.log,LastTest.log}`.
