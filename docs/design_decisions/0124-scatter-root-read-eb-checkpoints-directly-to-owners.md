# Decision 0124: scatter root-read EB checkpoints directly to owners

## Context

Root-only checkpoint writing removes rank-multiplied output buffers, but a
restart that reads or broadcasts the complete numerical hierarchy on every rank
would reintroduce the same memory window in the opposite direction.

## Decision

Read the established formatted multipatch checkpoint only on a caller-selected
root. Validate the reconstructed coarse geometry and patch topology against the
caller's replicated restart template. Copy entities owned by the root and send
one packed state/temperature payload for every root tile or child owned by a
different rank.

Require every non-root complete input array and root patch set to be empty at
the direct-scatter boundary. Broadcast only the small checkpoint clock metadata
and final read status. Publish the sparse candidate, metadata, and sender-local
traffic count only after collective validation succeeds.

## Consequences

Formatted restart no longer needs an all-rank numerical-field replica, and the
same serial checkpoint schema can target any valid current owner map. Traffic is
one root send per remotely owned entity.

The caller still provides the replicated geometry and full patch-set template.
Checkpoint-driven topology reconstruction without replicated child numerical
arrays requires a geometry-only descriptor and remains subsequent work.
