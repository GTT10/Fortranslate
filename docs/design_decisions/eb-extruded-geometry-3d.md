# Extrude validated 2D level-set sections into 3D EB metrics

`extrude_eb_geometry_3d(section,axis,layers,lower,upper,geometry,ok)` reuses the
existing 2D piecewise-linear level-set geometry instead of independently clipping
the same cross-section in 3D. The section must pass its existing validator.
Axis 1 maps section coordinates to (y,z), axis 2 to (z,x), and axis 3 to (x,y).
This cyclic convention preserves the solid-to-fluid normal orientation.

Volume fractions and transverse centroids repeat through the extrusion. Side-face
fractions and centroids copy section edges, while extrusion-normal faces copy the
section's fluid area and volume centroid. Boundary area and integrated normals
scale by layer thickness; boundary centroids use each layer's midpoint. Candidate
geometry is validated before publication, preserving the caller on rejected
inputs or invalid constructed metrics. Allocation uses the existing geometry
allocator; allocation exhaustion is not a recoverable status guarantee.

A cell may contain two differently oriented facets from the existing 2D triangle
split. Its scalar wall area is then larger than the magnitude of its integrated
normal. The 3D validator now uses the established 2D convention: a positive
integral magnitude no larger than area, parallel to the unit representative
normal, together with unchanged discrete face-closure checks. Tolerances are
unchanged. Planar geometry retains equality. The existing 3D pressure source
already uses integrated normals, so pressure closure requires no flux changes.

This supports prismatic and faceted cylindrical sections invariant along one
coordinate. It does not construct arbitrary 3D surfaces or resolve geometry
missing from the supplied 2D node samples. Existing 2D ambiguity/degeneracy checks
remain authoritative. CLI/restart geometry selection, general state redistribution,
curved-wall transport qualification, AMR and MPI EB integration are separate work.


## Outward-neighbor weighted StateRedist

The single-receiver weighted StateRedist no longer requires an axis-aligned wall.
For a small cell, it considers each nonzero solid-to-fluid normal component,
requiring an in-domain neighbor connected by a positive face aperture and with
fluid volume fraction at least the requested target. The largest eligible normal
component wins, with x/y/z tie order. This uses the physical unit normal, including
anisotropic cells. A receiver can be a regular cell or a sufficiently large cut
cell. It does not require its own merge at that target.

The original overlap counts, self/neighbor weights and component-wise conservative
averaging remain unchanged. An axis plane therefore retains its original receiver
and arithmetic. If no eligible single neighbor exists, the existing failure result
is retained; no chain, multi-neighbor region or traversal through solid is added.
This is a bounded extension, not general small-cell agglomeration. The coupled
library driver now has a faceted-cylinder test through hydro, transport and reaction;
application geometry/restart configuration is still separate work.


## Application configuration and restart

The common EB3D configuration now selects plane (default) or extruded geometry.
The latter reads a transverse nodal level-set file, validates its mesh dimensions
and finite samples, builds a private 2D section, then publishes validated extruded
geometry. The section path must differ lexically from output/checkpoint/restart
paths. Missing values are initialized to quiet NaN before list-directed input,
so incomplete/null records cannot become unintended geometry.

The existing checkpoint axis record stays x/y/z for planes and carries
`extruded:x`, `extruded:y` or `extruded:z` for extrusion. This prevents mode/axis
aliasing without changing the planar record layout. Existing metric-by-metric
geometry comparison binds restart to the section's constructed geometry, not its
file path. Standalone continuation still requires the section input. Geometry
matching tolerances are unchanged; raw file hashing is not introduced.

For each momentum direction the application checks whether any embedded normal
integral is nonzero. A wall can exert force in such directions, so they are
reported as momentum changes and excluded from the invariant test, just as the
single normal component was for an axis plane. Mass, energy, elemental and other
existing gates retain their tolerances. No integrated wall-force ledger or
curved-wall physical validation is implied.
