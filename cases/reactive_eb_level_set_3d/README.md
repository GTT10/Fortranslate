# Nodal 3D embedded-boundary application

This small stationary spherical-solid fixture runs fixed/selected chemistry,
molecular transport and StateRedist hydro with a fully three-dimensional nodal
level set. It checks integration and restart, not a calibrated sphere-flow problem.

From the repository root:

```sh
cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_3d
cd cases/reactive_eb_level_set_3d
../../build/selected-tree-mpi-debug/pelef_reactive_eb_3d reference.nml
```

The existing `reactive_eb_3d` namelist accepts:

```fortran
geometry_type='level_set'
level_set_file='level_set.dat'
```

The first record gives the three **cell counts** nx, ny, nz, matching the mesh.
The remaining values describe (nx+1)*(ny+1)*(nz+1) nodes with x varying fastest,
then y, then z. All values must be finite; positive denotes fluid. Coordinates
come from the configured Cartesian bounds and grid. File paths are relative to
the process working directory. The supplied 10x8x6 mesh has 11x9x7 node samples.
Trailing nonblank records are rejected. Plane/extrusion controls are unused in
this mode; keep their defaults when restarting.

Each cell is split into six tetrahedra with linear interpolation. Shared face
triangulations agree. This produces a faceted approximation, not an exact curved
surface. Cells with disconnected fluid regions, zero/degenerate geometry or
invalid face closure reject. Features absent from the nodal samples cannot be
reconstructed. StateRedist still requires a sufficiently large outward neighbor;
more complex small-cell clusters may reject during stepping.

Use `checkpoint_file`, `checkpoint_interval_steps`, `stop_after_checkpoint` and
`restart_file` as in the planar/extruded frontend. Geometry identity is `level_set`
in the checkpoint's geometry selector. Continuation rereads the requested nodal
file and compares its full constructed geometry against checkpoint metrics using
existing tolerances. The file may be relocated; a changed shape rejects before
output replacement. Raw nodal file bytes are not checkpoint identity.

Momentum components subject to wall forces are reported as changes; existing
mass/energy/element gates remain. This is not an integrated force ledger. The
selected frontend uses the same geometry interface and named composition inputs.
MPI/AMR EB integration, larger-neighborhood stabilization, convergence and physical
validation remain separate work.
