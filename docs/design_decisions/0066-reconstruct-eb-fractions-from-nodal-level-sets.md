# Decision 0066: reconstruct EB fractions from nodal level sets

## Context

Embedded-boundary flow needs a geometry representation before cut-cell fluxes,
small-cell stabilization, or wall conditions can be qualified. Coupling a new
geometry algorithm directly into the reactive CTU update would make errors in
volume fractions, face apertures, and flow physics difficult to isolate.

## Decision

Add a standalone two-dimensional Cartesian geometry module. A caller provides
finite level-set values on every mesh node; positive values define fluid. Split
each cell along the lower-left to upper-right diagonal, treat the level set as
affine on each triangle, clip its positive polygon, and sum normalized polygon
areas for the cell-volume fraction.

Compute each Cartesian face open fraction once from linear interpolation of
its two endpoint values. Classify a cell as covered or regular only when its
volume fraction lies within a fixed roundoff tolerance of zero or one;
otherwise classify it as cut. Keep all fractions and classifications in a
self-validating geometry value independent of hydro state.

## Consequences

Planar interfaces are represented exactly by the piecewise-affine geometry,
and smooth curved interfaces converge under node refinement. Neighboring cells
share the same face-fraction storage, preventing two reconstructions of one
aperture.

The fixed triangle diagonal is an intentional first geometry contract, not an
AMReX EB2 reproduction. The milestone does not yet compute cut-face centroids
or normals and does not change a fluid update. Flux divergence,
redistribution/small-cell stabilization, embedded-wall states, AMR, and MPI
remain later EB work.
