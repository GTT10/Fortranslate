# Decision 0130: select sparse EB timestep on tile owners

## Context

The sparse EB timestep selector gathers all root tiles to the first tile's
owner before evaluating hyperbolic and parabolic stability limits. Both limits
are cell-local reductions: they depend on the conserved state, temperature,
cell spacing, and transport data, but not on neighboring tile values. The
gather therefore transfers and allocates a complete root field without a
numerical dependency that requires it.

## Decision

Extract the exact EB geometry row band for every locally owned root tile and
evaluate its reactive CFL and optional molecular-transport timestep directly
on that tile's sparse state. Continue evaluating fine children only on their
owners and scale their limits by the refinement ratio. Skip fully covered
entities, then select the global coarse interval with one communicator minimum.

Retain the public transfer counter for compatibility but require it to remain
zero. Reject an entirely inactive hierarchy instead of accepting the finite
`huge()` sentinel as a timestep. Preserve collective validation and publish dt
and traffic only after every owner succeeds.

## Consequences

Every accepted-step stability selection in the sparse public clock now has no
root numerical-field communication or complete-root allocation. Rank-count
parity is exact because the global operation remains a minimum over the same
cell-local candidate limits.

SSPRK2 transport flux construction and StateRedist do require neighboring
values and still use selected-root gathers. Their decomposition needs explicit
halo exchange and remains separate work.
