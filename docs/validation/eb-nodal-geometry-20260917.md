# Nodal 3D embedded-boundary geometry — 2026-09-17

## Source and scope

Base `605ec22c17ce8b6cd83a4ed81d280e8713754623`, local `main`.
Adds triangle/tetrahedron clipping and a transactional 3D nodal level-set geometry
builder. Existing pressure fluxes, geometry tolerances, mechanisms and generated
sources are unchanged. Fixed and selected EB3D runtimes include the clipping module.

SHA-256:

- `src/eb/eb_polyhedron_clip_mod.F90`:
  `88ce5c52b4c50b1a7726cf551ab53d8631a4be905a1ab4d870fb74291955747b`
- `src/eb/eb_geometry_3d_mod.F90`:
  `2689eac3ef56abfd9d3cf4e2c4ffe6106d939c6d1d4773d2e6261b151503fd01`
- `tests/unit/test_eb_geometry_3d.F90`:
  `cb81df0795acc5a6447752207a1a307417f3322d81bfe20dff891aa0c1127142`
- `tests/unit/test_eb_reactive_wall_flux_3d.F90`:
  `7082e30b125ac61ad46ea51e19499fd948d26a5cc8998ff85a5e5e708755f676`

## Commands and completed results

GNU Fortran 13.3.0 Debug; existing `build/selected-tree-mpi-debug` configured with
OpenMPI 4.1.6. These geometry/pressure tests execute serially.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_geometry_3d test_eb_reactive_wall_flux_3d pelef_selected_reactive_eb_3d_runtime -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_(geometry|reactive_wall_flux)_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed core and selected runtime builds passed. Final CTest **2/2 passed**, no
failures/skips, 0.10 s total (geometry 0.02 s, wall flux 0.05 s). Earlier focused
passes preceded addition of pressure and disconnected-cell checks; the final pass
includes both. Provenance, contract and whitespace checks passed. No Python source
changed. Existing compiler warning categories remain.

## Evidence and limits

For x+y+z>2.13 in a unit cube, tests compare the analytic corner-tetrahedron volume,
centroid and oblique interface area on a 10x8x6 grid. A sampled spherical boundary
and its sign-reversed complement check per-cell volume partition, equal interface
area, opposite integrated normals and discrete face closure. Nonfinite inputs and
disconnected positive-volume cells preserve prior geometry on rejection. Existing
plane/extrusion tests remain. The real 3D pressure-source/divergence kernel balances
uniform Cartesian pressure against the fully 3D faceted sphere at the existing
cell-integrated tolerance.

The sphere is a tetrahedral piecewise-linear reconstruction, not an exact sphere.
Zero/degenerate metrics and multi-component fluid cells are unsupported. No nodal
3D CLI/restart run, curved 3D coupled evolution, convergence, EB-AMR-MPI integration,
full matrix, installed audit, external solver parity or physical validation was run.
Logs: `build/eb-nodal-geometry-checks/{build,test-build,flux-build,final-build,ctest}.log`.
