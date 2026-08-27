# Decision 0159: Drive the EB tree clock on owners

## Context

The sparse arbitrary-depth tree has owner-local timestep selection and an
atomic owner-local full-physics step. A public clock must combine them without
allowing rank-dependent time controls, an uncommitted step, or diagnostics that
advance farther than the accepted numerical state.

## Decision

Establish communicator agreement on time, target, step counts, step ceiling,
CFL values, physics controls, and data extents before entering the loop. Before
every attempted step, recompute the sparse owner-local stability limit and clip
it to the remaining target interval. Advance a private sparse full-physics
candidate and commit fields, time, step, minima, advances, timestep-node count,
and transfer categories together.

Keep already committed steps when a later timestep or physics operation fails
or the step ceiling is reached. Publish the requested target time exactly only
on successful termination.

## Consequences

Normal time advancement never materializes the arbitrary-depth numerical tree.
Clock diagnostics identify exactly the committed owner work and traffic. An
initial rank mismatch or zero-step ceiling mutates nothing, while a later
failure preserves a consistent accepted prefix. Dynamic tagging and
checkpoint/restart can be added around this commit boundary separately.
