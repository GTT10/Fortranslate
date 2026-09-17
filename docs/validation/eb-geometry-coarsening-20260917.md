# Conservative 3D EB geometry coarsening — 2026-09-17

## Implementation and boundary

Source: `009f1812ae7c74b1e130b2dfec218a2868e57b54`, based on
`37b2425c992145f1b261a7fea64be9b9842ce805` on local `main`.

`coarsen_eb_geometry_3d(fine, ratio, coarse, ok)` constructs a full-domain
coarse geometry from an integer subdivision of an existing fine geometry.
It preserves fluid volume and its first moments, sums wall area, wall first
moments and integrated normals, and aggregates only the coincident exterior
fine faces for coarse face apertures and tangential centroids. Internal faces
cancel in the discrete normal closure. Physical bounds are unchanged; cell and
face centroid offsets are transformed to coarse normalized coordinates.
Ratio one copies the geometry exactly.

Each coarse cell must contain one connected fluid region through positive fine
face apertures. Disconnected regions and wall-normal cancellation that prevents
a representative unit normal are rejected. The existing geometry validator and
its tolerances are unchanged. Invalid ratios, incompatible dimensions, invalid
input geometry or unrepresentable coarse geometry leave the destination intact.
Allocation exhaustion follows the existing geometry-builder behavior and is not
a recoverable failure contract.

This library operation is a prerequisite for consistent EB coarse/fine metrics.
It does not yet connect EB to patch-tree stepping, state restriction, reflux,
regridding or MPI. It does not resample a level set or smooth away unresolved
closed surfaces. General EB-AMR integration remains open.

## Independently completed checks

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; matching existing build
`build/selected-tree-mpi-debug`, MPI enabled and spray disabled.
Fixed core and selected EB runtime compiled. No MPI execution or application
matrix was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_geometry_3d test_eb_reactive_wall_flux_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_geometry_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_reactive_wall_flux_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_geometry_3d|unit_eb_reactive_wall_flux_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 PASS, 0 failures, 0 skips, 0.15 s**:

- Geometry: 0.06 s. Ratio-two and ratio-three restriction of a 12-cubed oblique
  plane in an anisotropic translated box matches independent coarse construction
  for all volume, face, centroid and wall metrics at the existing tolerance.
  A faceted curved surface preserves volume, area and wall first moments.
  Invalid input/ratio, disconnected fluid regions and a closed surface with
  cancelling normals reject transactionally. Ratio-one identity is checked.
- Wall flux: 0.06 s. The existing uniform-pressure divergence check additionally
  uses a 20-by-16-by-12 nodal sphere restricted to 10-by-8-by-6, preserving
  pressure balance at the existing tolerance.

Provenance and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
A public-declaration indentation fix followed the test run; numerical source was
unchanged. Compiler logs include existing warning categories.

Logs: `build/eb-coarsen-checks/{build,test-build,final-build,ctest}.log`.
The earlier geometry-only pass (`geometry.log`) is not added to the final count.
This is focused numerical evidence, not external solver parity, convergence,
physical validation or overall project completion.
