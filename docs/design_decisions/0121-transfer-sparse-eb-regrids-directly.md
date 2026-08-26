# Decision 0121: transfer sparse EB regrids directly

## Context

The sparse MPI EB clock and tag planner are owner-local, but their final
two-level topology transaction still materializes complete root and child
fields on every rank. That correctness bridge preserves serial regrid semantics
at the cost of rank-multiplied storage and numerical-field broadcasts during
each dynamic topology change.

## Decision

Preserve the serial operation order without constructing a replicated field
tree. First apply the existing targeted old-fine-to-root average-down to a
private sparse candidate. Build the new replicated geometry/template and its
deterministic owner map. Send each averaged root tile once to every distinct new
child owner that does not already own it; those owners perform PCM
prolongation locally.

Enumerate same-ratio old/new overlap rectangles from replicated patch metadata.
Validate cell types and volume fractions collectively, then copy a rectangle
locally when its owner is unchanged or send it once from the old owner to the
new owner. Recover active-cell temperatures on new child owners. Publish the
new distribution, sparse fields, geometry template, and restriction,
prolongation, and overlap transfer counts only after every rank accepts the
complete candidate.

## Consequences

Explicit, tagged, and scheduled two-level regrids remain globally single-copy
for numerical fields and retain the established serial average-down, PCM,
overlap, and temperature-recovery behavior. Failure leaves the original
distribution and payloads untouched and publishes zero traffic for an
uncommitted transaction.

Geometry and compact topology descriptors remain replicated intentionally.
The root may be assembled on more than one distinct new child owner, so root
traffic grows with participating owner count rather than child count. Sparse
checkpoint/output, arbitrary-depth dynamic EB topology, and decomposition of
the level-wide root physics kernel remain separate work.
