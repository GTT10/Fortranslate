# Extruded 3D embedded-boundary geometry — 2026-09-17

## Source

Base `4bd974d3bf3333b487ae75bd9ced72bf2e3bb79d`, local `main`.
Adds transactional three-axis extrusion of valid 2D level-set sections and the
faceted integrated-normal geometry contract. Pressure flux implementation,
mechanisms and numerical tolerances are unchanged. Selected EB3D runtime now
includes the existing 2D geometry dependency.

SHA-256:

- `src/eb/eb_geometry_3d_mod.F90`:
  `c8ba7d8005369d43f1975cfc849c2fa49a11ef96d046f9668b93dbeb3febf0a8`
- `tests/unit/test_eb_geometry_3d.F90`:
  `6c7b283e25ecfd562cf379c20a5061ed5210a222898c3705ef74865160bd2da9`
- `tests/unit/test_eb_reactive_wall_flux_3d.F90`:
  `5836646490a86f3b47cad2c33394e525a09ea14d9bbfb8dc3788216235344111`

## Commands and completed checks

GNU Fortran 13.3.0 Debug in existing `build/selected-tree-mpi-debug` (OpenMPI 4.1.6
configuration; these focused geometry/flux tests are serial).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_geometry_3d test_eb_reactive_wall_flux_3d pelef_selected_reactive_eb_3d_runtime -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_geometry_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_(geometry|reactive_wall_flux)_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed core and selected-runtime builds passed. Final CTest **2/2 passed**, no
failures/skips, 0.05 s total (0.01 s each). An initial focused pass preceded the
addition of analytic centroid assertions; the recorded final pass includes them.
Provenance, contract and whitespace checks passed. No Python source changed.
Logs: `build/eb-extruded-geometry-checks/{build,centroid-build,ctest}.log`.

All three extrusion axes check the analytic volume, centroid and wall area of
the triangular fluid section x+y>1.13 in a unit square, extruded over length 2.5.
A sampled circular solid checks faceted section-to-3D volume/area and discrete
normal closure. Invalid axis/layer/bounds/section inputs preserve prior geometry;
area smaller than the integrated normal is rejected. Existing planar corruption
and alignment tests remain. A faceted cylindrical wall uses the real pressure
source/divergence with uniform Cartesian pressure fluxes and satisfies the
existing cell-integrated pressure-balance tolerance.

## Limits

Library geometry and pressure balance only; no nonplanar application run,
curved-wall transport/state-redistribution qualification, AMR/MPI geometry,
full matrix, installed audit, external solver parity or physical validation.
The circle is the 2D piecewise-linear approximation, not an exact curved surface.
