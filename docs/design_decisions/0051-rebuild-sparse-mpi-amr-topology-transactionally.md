# Decision 0051: rebuild sparse MPI AMR topology transactionally

## Context

Sparse patch payloads can migrate between owners while a hierarchy is fixed,
and all physics operators now run directly on that storage. A changed regrid
plan also changes patch geometry, counts, parent relationships, and the
deterministic load-balanced owner map.

## Decision

Provide one explicit-plan sparse regrid transaction. Materialize the current
sparse hierarchy collectively, apply the qualified serial conservative rebuild,
and require every rank to report the same change decision and overlap-transfer
count. If topology changes, initialize a new owner distribution and scatter the
rebuilt tree to exactly those owners.

Treat the sparse solution and owner distribution as one commit unit. An
unchanged plan updates only evaluation bookkeeping. Any invalid plan,
inconsistent result, distribution failure, or scatter failure restores the old
sparse solution and returns the old distribution.

## Consequences

Explicit patch-tree topology changes preserve serial overlap, conservation,
bookkeeping, and rollback semantics while returning to globally single-copy
persistent storage. Patch additions can alter the deterministic owner map, so
new owners receive the rebuilt authoritative payloads.

The transition temporarily materializes a full reactive replica on every rank.
Direct tag-driven planning, owner-to-owner overlap transfer, and point-to-point
communication remain separate work.
