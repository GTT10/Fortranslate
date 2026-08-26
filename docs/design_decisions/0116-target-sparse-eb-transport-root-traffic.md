# Decision 0116: target sparse EB transport root traffic

## Context

Direct sparse SSPRK2 transport keeps fine children owner-local, but each Euler
stage materializes the complete root on every rank. It broadcasts the updated
root state, temperature, diffusive fluxes, and every cumulative reflux
correction. The final SSPRK2 blend also materializes both root candidates and
broadcasts the blended field. Unrelated ranks neither advance root physics nor
own a fine child or the affected root row bands.

## Decision

For each Euler stage, gather nonlocal root tiles only to the root physics owner.
Advance the existing root transport and StateRedist kernels there, then send one
packed start/update/temperature/flux bundle to each distinct remote child
owner. Preserve child ordering by sending the current corrected root to a
remote child before its subcycles and returning the refluxed root afterward.
Scatter only the final row bands to remote root tile owners.

For the final SSPRK2 blend, gather the start and second-Euler root candidates
only to the physics owner, recover the blended temperature there, and scatter
only the final row bands. Compute the EB boundary-change vector on the physics
owner, broadcast only that `nvar` vector, and apply its conservation correction
directly to locally owned root tiles.

Count each packed point-to-point MPI message as one root payload transfer.
Publish transfer counts, Euler counts, limiter minimum, and sparse state only
after the complete SSPRK2 transaction succeeds.

## Consequences

Unrelated ranks no longer allocate or receive complete transport root fields.
Communication scales with nonlocal root tiles, distinct child owners, and
remote children while retaining root/child update ordering, serial parity, and
transactional rollback.

The root transport stencil and StateRedist remain one level-wide operation on
the root physics owner. Root physics decomposition, public sparse time-loop
control, dynamic sparse topology, checkpointing, and output remain future work.
