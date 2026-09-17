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
