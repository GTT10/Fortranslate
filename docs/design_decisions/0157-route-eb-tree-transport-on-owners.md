# Decision 0157: Route EB tree transport on owners

## Context

Recursive SSPRK2 molecular transport needs two complete Euler traversals, each
with parent-time exterior context, fine-flux accumulation, reflux, ordered
average-down, StateRedist, and subtree conservation closure. Materializing the
tree between stages would erase the sparse-ownership guarantee established for
hydro.

## Decision

Execute each recursive transport node only on its sparse owner. Reuse the
qualified direct hydro route for compact parent context, child fine-flux
return, parent-owner registers, reflux state round trips, and ordered
average-down. Keep transport flux/RHS evaluation and StateRedist on the owner.

Retain the accepted tree privately through both Euler stages. Blend the start
and second-stage fields and recover temperature locally on every owner, then
perform one final direct deepest-first restriction. Reduce the limiter minimum
across the communicator and publish fields, advances, transfers, and limiter
only after collective validation.

## Consequences

No nonowner allocates a numerical node during recursive transport. Two Euler
traversals cost twice the hydro-style `r + 4` remote-edge schedule, and the
post-blend hierarchy synchronization adds one payload per distinct-owner
relation. Shared-owner work remains local. The result follows the serial
SSPRK2, reflux, conservation, and restriction order within qualified
roundoff. Full-physics composition and its public clock remain separate.
