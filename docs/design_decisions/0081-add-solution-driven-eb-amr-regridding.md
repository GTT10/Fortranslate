# Decision 0081: add solution-driven EB AMR regridding

## Context

The runnable reactive EB AMR application owned a qualified static fine
rectangle, but could not move that rectangle when the resolved temperature
structure lay elsewhere. Replacing a patch must not lose its composite
conserved state or overwrite same-resolution data that remain refined.

## Decision

Tag active, strictly internal root cells from the largest four-neighbor
temperature jump. Require both an absolute jump floor and a normalized relative
threshold. Enclose all tags in one buffered, minimum-size rectangle that stays
one coarse cell away from every physical boundary.

Treat topology replacement as one transaction. First average the complete old
fine patch into the root. Piecewise-constant prolongation then initializes the
new rectangle from that synchronized root. Finally, copy the exact old fine
state and temperature into every overlapping global fine index and recover the
EOS on all active new fine cells before committing either level.

Run this planner at initialization and after every configured number of accepted
coarse steps. Keep the existing fine patch when no cells are tagged or when the
planned bounds are unchanged. Report the number of committed regrids.

## Consequences

The serial two-level EB application can now move and resize one rectangular
fine patch around a temperature feature without changing the composite
conserved integral. Retained overlap is not diffused by repeated prolongation.

The policy intentionally does not remove the only fine patch when tags vanish.
It also does not create multiple patches, deeper levels, boundary-touching
patches, EB AMR chemistry or transport, checkpoints, or MPI ownership.
