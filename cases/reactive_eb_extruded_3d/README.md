# Extruded 3D embedded-boundary application

The fixed and configure-time selected EB3D frontends accept an extruded nodal
level-set section. The supplied small case uses a faceted circular solid,
StateRedist, chemistry and molecular transport. It is an integration fixture,
not a calibrated cylinder-flow experiment.

Build and run from the repository root:

```sh
cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_3d
cd cases/reactive_eb_extruded_3d
../../build/selected-tree-mpi-debug/pelef_reactive_eb_3d reference.nml
```

The existing `reactive_eb_3d` group adds:

```fortran
geometry_type='extruded'
extrusion_axis='z'
section_file='section.dat'
```

The default remains `geometry_type='plane'`. Extrusion uses cyclic transverse
coordinates: x -> (y,z), y -> (z,x), z -> (x,y). Domain bounds and grid extents
come from the ordinary 3D inputs. `plane_axis` and `plane_position` do not build
extruded geometry; retain their defaults across a restart.

`section_file` is relative to the process working directory. Its first record
contains the two transverse **cell counts**, matching the 3D mesh. Remaining
records contain (first count + 1) * (second count + 1) finite nodal values, with
the first coordinate varying fastest. Positive values denote fluid. The supplied
file has an 11-by-9 node grid for a 10-by-8 cross-section. The existing 2D
piecewise-linear builder defines the faceted interface and rejects unsupported
ambiguous/degenerate sections. Trailing nonblank records are rejected.

Use the existing `checkpoint_file`, `checkpoint_interval_steps`,
`stop_after_checkpoint` and `restart_file` controls. Restart needs a section file
that reproduces the saved geometry. It can be relocated; the path is not physical
identity. Geometry mode/axis and all stored geometry metrics are checked using
existing checkpoint tolerances. A changed section that changes these metrics is
rejected before output replacement. Old plane checkpoint records retain their
axis identity. Older executables do not support extruded checkpoint identities.

Momentum components on which the embedded wall can exert pressure forces are
reported as changes, rather than treated as invariants; other existing conservation
gates remain. This is not a wall-impulse ledger. Single-neighbor StateRedist still
requires an outward open-face neighbor that meets its volume target; narrow
features requiring multi-cell agglomeration can fail. The geometry is invariant
along the extrusion direction. Fully 3D surfaces, MPI/AMR EB, moving walls and
physical validation remain separate work.
