# Decision 0056: send parent states only to child owners

## Context

Final sparse ghost refresh uses each parent state to reconstruct ghosts for its
children. The previous correctness schedule broadcast the complete parent state
to every rank even though only ranks owning those children consumed it. Several
children may share a remote owner and can reuse the same parent state.

## Decision

For each parent, derive the distinct ranks that own at least one child. Keep the
authoritative state local on the parent owner and send one complete state to
each distinct remote child owner in rank order. Each recipient reuses that copy
for all of its local children. Ranks owning neither the parent nor a child skip
allocation and communication.

Count a successful send on the parent owner for every remote recipient, and
compare the communicator sum with an independent traversal of the owner map.
Retain the surrounding collective acceptance boundaries and transactional
rollback.

## Consequences

Final parent-to-child ghost refresh no longer scales traffic with communicator
size or duplicate children on a recipient. It still sends the complete parent
state; narrowing that payload is separate optimization work.

Recursive hydro/transport interval state and face-flux distribution, level
counter synchronization, flux reconciliation, and topology-changing overlap
transfer remain collective.
