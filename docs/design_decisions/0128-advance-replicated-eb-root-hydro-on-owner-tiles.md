# Decision 0128: advance replicated EB root hydro on owner tiles

## Context

The replicated MPI-owner hydro transaction assigns root storage by y-tile but
still chooses the first tile's owner to advance the complete root level. That
rank computes every root face, divergence, and StateRedist update, then
broadcasts the complete state, temperature, x-flux, and y-flux arrays. Fine
children run on their owners, so the remaining level-wide root kernel is the
hydrodynamic compute bottleneck.

## Decision

Advance one bounded y-band per root tile on that tile's owner. Extend the band
by at most six cells above and below the owned rows. Extract all cell, face, and
embedded-boundary metrics with local array bounds and invoke the established
reactive EB level kernel. Six guard rows conservatively cover the combined
reconstruction, face-centroid interpolation, overlapping StateRedist
neighborhood, and second-order slope dependencies.

Publish state, temperature, and x-face flux only for the owned cell rows.
Assign every y-face to the tile immediately above it, except that the final
tile also owns the upper physical boundary. Sum zero-filled rank contributions
to assemble complete replicated outputs. Count local tile advances and actual
halo-band cells, but publish those counters only if the entire root, child,
reflux, and average-down transaction commits.

## Consequences

Root hydrodynamic work is distributed across the existing root owners, and no
rank advances the complete root solely to broadcast four result arrays. The
same serial kernel and bounded-overlap construction preserve face and
StateRedist parity without treating tile edges as physical boundaries.

This decision applies to the replicated owner path. Its complete root inputs
and assembled outputs remain replicated for compatibility with existing child
exterior construction and deterministic reflux. The sparse path still gathers
root tiles to one physics owner; direct sparse halo exchange and sparse result
routing remain subsequent work.
