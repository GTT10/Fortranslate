# Decision 0069: preserve the integrated EB normal

## Context

A cut cell can contain one affine interface segment from each of its two
triangles. Averaging their normals by length and renormalizing is useful for a
representative unit wall normal, but multiplying that unit vector by total
interface length does not recover the sum of the two oriented segment vectors
when the piecewise-linear interface bends inside the cell. That discrepancy
breaks constant-pressure balance on curved embedded boundaries.

## Decision

Store both the representative solid-to-fluid unit normal and the integrated
normal `integral(n ds)`. Sum each segment's unit normal times its physical
length before normalization. Use the unit normal for local wall orientation and
the integrated normal for the cell-integrated pressure force.

Form the Cartesian part of the cut-cell divergence from shared face fluxes
multiplied by their open fractions and physical face lengths. Add the wall
source and divide the complete balance by `kappa*dx*dy`. Covered cells have
exactly zero right-hand side, and rejected inputs leave the full output zero.

## Consequences

Constant pressure now satisfies the discrete surface-vector identity for both
planar and piecewise-linear curved cells. This refines Decision 0068: the
`p*n*A` expression is retained for a single planar normal, while aggregated
segments use `p*integral(n ds)`. Small-cell stability and flux construction
beside covered cells are still unresolved.
