# Decision 0182: Regrid a fixed three-level parent transactionally

## Context

The fixed three-level driver could move its middle-to-finest rectangle, but
its root-to-middle rectangle remained fixed. Moving that parent invalidates
the finest patch coordinate system, so publishing the middle regrid before a
replacement finest patch exists would expose an incomplete hierarchy. The old
finest state also owns the composite solution over part of the middle level.

## Decision

Add one library transaction for a root-tagged parent regrid. Restrict the old
finest state into the old middle first, then use the established two-level
regrid to average the synchronized middle into the root, prolong the new
middle, and retain its same-resolution overlap. Plan a new finest rectangle
inside the required two-cell middle margin and initialize it through the
configured conservative prolongation method.

Compute the three-level composite conserved integral before and after the
candidate rebuild. Commit the root and both refined state/temperature arrays,
geometries, and patch descriptors together only when the topology, EOS
recovery, and conservation check all succeed. Missing interior finest tags or
any transfer failure rejects the complete candidate.

## Consequences

The fixed-depth library can now relocate or resize a parent patch without
losing the old finest contribution or exposing a two-level intermediate
state. Same-resolution middle overlap is retained. Old finest subcell detail
is represented through its conservative restriction before the replacement
finest patch is prolonged. Scheduling this operation in the public time loop
and extending the fixed-depth checkpoint topology remain separate lifecycle
work; arbitrary-depth serial and sparse-MPI trees already regrid every parent.
