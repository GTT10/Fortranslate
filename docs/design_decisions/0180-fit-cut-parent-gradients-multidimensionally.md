# Decision 0180: Fit cut-parent gradients multidimensionally

## Context

The first conservative cut-parent prolongation measured independent x and y
secants. That preserves variation next to a simple embedded boundary, but a
fluid-volume centroid can shift in both coordinates. Treating the difference
to an axial neighbor as a pure coordinate derivative therefore mixes the
normal and tangent gradients and does not exactly reproduce a general affine
field at cut parents.

## Decision

Fit one component-wise linear gradient from the active coarse fluid-volume
centroids in the surrounding 3-by-3 cells. Axial samples require their shared
face to be open. Diagonal samples require at least one open two-face path
through an active intermediate cell.

Solve the two-by-two normal equations when they are full rank. For a rank-one
stencil, use the minimum-norm gradient in the resolved direction; use zero
gradient when the stencil has no direction. Limit predictions at connected
coarse centroids to the active local component envelope, then retain the
existing fine-child envelope limiter.

## Consequences

Affine conserved fields evaluated at fluid-volume centroids reproduce exactly
at active fine centroids when the child predictions remain inside the local
envelope; the qualified diagonal-plane case exercises the interface tangent.
The volume-weighted zero-mean correction keeps cut-parent average-down
conservative. Disconnected cells do not influence the fit. The limiter may
reduce nonlinear or outward-normal gradients. This is not an exact AMReX
interpolation or a quadratic cut-cell reconstruction.
