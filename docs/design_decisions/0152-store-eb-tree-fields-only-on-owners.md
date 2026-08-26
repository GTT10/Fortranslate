# Decision 0152: Store arbitrary-depth EB tree fields only on owners

## Context

The arbitrary-depth MPI tree has deterministic node ownership and exact rank
accounting, but its first publication operation retains complete numerical
fields on every rank. Owner-local timestep and recursive physics would not be
meaningfully sparse while nonowners still allocate every node.

## Decision

Keep topology and the distribution descriptor replicated. Introduce a sparse
numerical tree whose state and temperature arrays exist exactly on the rank
that owns each node. Treat reconstruction of a complete numerical tree as an
explicit materialization operation.

For a distribution change, allocate a private sparse candidate from the new
map. Copy unchanged-owner nodes locally and transfer changed nodes directly
from old owner to new owner with ordered state and temperature messages. Check
that every rank holds identical old and new distribution descriptors before
communication, validate the complete sparse candidate collectively, and only
then replace the accepted tree.

## Consequences

Nonowners no longer retain persistent numerical fields. Replicated consumers
remain available through an explicit, measurable boundary, while normal
ownership changes avoid an all-rank field broadcast. Invalid or inconsistent
owner metadata rejects before transfer and preserves the accepted sparse tree.

Stable-step selection and recursive physics still consume the serial complete
tree and therefore remain the next owner-local routing boundaries.
