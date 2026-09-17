# Build 3D EB metrics from nodal level sets

`build_level_set_eb_geometry_3d(phi,lower,upper,geometry,ok)` constructs cell and
face metrics from a three-dimensional array of nodal values. Positive values are
fluid. Each Cartesian cell uses the same six-tetrahedron decomposition around
its lower-to-upper diagonal. Shared faces therefore use matching triangle
partitions. Interpolation is linear within each tetrahedron; a curved implicit
surface becomes a faceted approximation.

`eb_polyhedron_clip_mod` clips triangles and tetrahedra against zero level set.
Fluid face polygons and the ordered interface cap form a convex polyhedron.
Triangular pyramids from an interior point provide volume and first moments;
interface facets provide scalar area, area moment and solid-to-fluid normal
integrals. Input orientation does not affect pyramid volumes. Computation uses
cell-local physical coordinates to avoid subtraction of large global moments.
Nodal values are scaled before interpolation. Complete-cell metrics aggregate
these extensive quantities; Cartesian face fractions/centroids come from the
matching clipped triangles. Existing facet-normal and closure validation remains
unchanged.

Positive vertices connected within each tetrahedron define the fluid connectivity
graph. A cell with disconnected positive-volume fluid components is rejected,
because one state per cell cannot represent them independently. All-zero triangles
or tetrahedra, unresolved zero-normal interfaces, and metrics inconsistent with
existing regular/covered-cell classification reject. This does not infer subcell
features absent from the nodal samples. No tolerance is loosened to admit a mesh.

The full grid is constructed privately and published only after validation.
Invalid inputs or unsupported geometry preserve the prior output. Existing
allocation behavior remains: memory exhaustion is not a recoverable status
contract. This increment exposes a library geometry builder and pressure-balance
checks. Nodal-file application/restart wiring, coupled 3D curved hydro/transport,
small-cell stabilization across larger neighborhoods and EB-AMR-MPI integration
remain separate work.


## Application and restart integration

`geometry_type='level_set'` selects a three-dimensional nodal file through
`level_set_file`. Dimensions must match the configured grid; missing values are
initialized to quiet NaN before input so incomplete records cannot become geometry.
The builder publishes validated metrics before initialization or checkpoint loading.
Input geometry paths must differ lexically from output/checkpoint/restart paths.

The existing geometry-selector checkpoint record stores `level_set`, distinct from
plane axes and extruded selectors. Existing full metric comparison binds restart
to the requested reconstructed geometry and retains its tolerances. The source file
may move; raw file bytes are not identity. No new checkpoint layout is needed.
The common fixed/selected driver then uses the already connected nonplanar
StateRedist, molecular transport, reaction and wall-momentum diagnostics.
