# Decision 0067: aggregate EB interface segments by physical length

## Context

Cell and Cartesian-face fractions identify cut cells but do not locate the
embedded wall or provide an orientation for a wall flux. The fixed-diagonal
piecewise-affine representation may produce one zero-contour segment in each
of the cell's two triangles.

## Decision

Extract the `phi=0` segment and normalized `grad(phi)` from each affine
triangle. Combine segment centroids and normals using physical segment length,
then normalize the accumulated normal so it points toward positive `phi`, the
fluid side. If both triangles report the same diagonal segment, retain it once.
Store zero metrics for regular and covered cells and validate cut centroids,
positive lengths, and unit normals with the rest of the geometry.

## Consequences

Planar wall location and orientation are exact, while curved-wall perimeter and
normal errors converge with refinement. The geometry can now support a cut-wall
flux, but this decision does not define that flux or a small-cell update.
