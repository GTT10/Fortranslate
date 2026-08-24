# Decision 0058: regrid explicit sparse topology directly

## Context

The explicit-plan sparse regrid API returned globally single-copy storage, but
its transaction first gathered every owned field patch into a complete tree on
every rank, ran the serial rebuild everywhere, and scattered the result. This
made topology changes scale with the complete hierarchy on every process.

## Decision

Replicate only compact hierarchy and owner-map metadata. Construct the candidate
hierarchy collectively, allocate only new local patches, and traverse it from
root to leaves. Each parent owner performs conservative prolongation and sends
one interior-state payload only when a child has a different owner.

Enumerate same-resolution old/new patch overlaps from hierarchy geometry. Copy
local segments directly and send packed state plus temperature cells from the
old owner to the new owner. Finish with the existing distributed average-down
and ghost-refresh paths, then publish the new solution and owner map together.

Expose local counts for cross-owner prolongation and overlap messages. An
independent application traversal must reproduce their communicator sums.

## Consequences

Explicit-plan regrid no longer creates a complete old or new field replica.
Unchanged plans remain metadata-only, and invalid plans retain exact rollback.
The tag-driven API still materializes a tree to derive its plan; distributed
owner-local tagging is a separate next step.
