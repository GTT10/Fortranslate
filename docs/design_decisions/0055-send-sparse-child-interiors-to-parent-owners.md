# Decision 0055: send sparse child interiors to parent owners

## Context

Sparse chemistry average-down and recursive hydro/transport synchronization
need every child interior on exactly one parent owner. The previous correctness
schedule broadcast each complete child interior through the communicator even
though unrelated ranks never consumed it.

## Decision

Use one shared ordered transfer helper for both average-down and physics
synchronization. If child and parent owners differ, send the contiguous child
interior directly from the child owner to the parent owner. If ownership is the
same, copy locally. Do not allocate an interior receive array on unrelated
ranks.

Retain deterministic relation, parent, and child traversal so the blocking
send/receive pairs match without additional replicated schedule state. Count a
successful remote transfer on the child owner during chemistry average-down,
then compare the communicator sum with the owner-map-derived cross-owner child
count.

## Consequences

Child interior communication now reaches only its consumer and the same helper
covers chemistry, hyperbolic reflux synchronization, and diffusive reflux
synchronization. Existing collective acceptance and rollback boundaries remain
unchanged.

Parent interval states and fluxes, parent-to-child ghost fill, sibling flux
reconciliation, and topology-changing overlap transfer still use collective
communication.
