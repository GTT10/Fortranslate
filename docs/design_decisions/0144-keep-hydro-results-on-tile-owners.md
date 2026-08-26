# Decision 0144: Keep hydro results on tile owners

## Context

After `0.151.0`, sparse hydro fluxes remain on their root tile owners and route
directly to child owners. Stage-start and stage-end state and temperature still
move to a root physics owner, which extracts child context, merges ordered
reflux corrections, and scatters final rows. Transport already demonstrates
that patch-plus-two state/support and ordered corrections can route directly
between the same owners.

## Decision

Retain stage-start, uncorrected stage-end, and current corrected hydro state and
temperature on each root tile owner. Reuse the qualified direct state-support
route to assemble patch-plus-two arrays on each child owner and extract the
four-edge start/end interpolation context there. Reuse the qualified direct
correction route to return refluxed support to every intersecting tile owner
before the next child advances.

After all children, commit each corrected root tile locally. Remove the remote
tile-result message, complete root hydro state assembly, root-owner ordered
support merge, and final row scatter. Retain the finite-band halo exchange,
direct interface-flux route, child-local subcycling/reflux, and hierarchy-wide
average-down. Delete the superseded private complete-root bundle, compact root-
to-child context, correction, and scatter helpers and their message tags.

## Consequences

Sparse hydro has no post-compute root physics owner and allocates no complete
root state, temperature, or flux result. Point-to-point traffic contains only
finite-band halos and direct state, flux, and correction fragments. Child order
remains deterministic because every correction completes before the next child
support is assembled. Qualification must preserve exact serial fields,
conservation, owner work, rollback, and message counts in Debug and Release at
one, two, four, and eight ranks. This decision does not claim a measured
speedup.
