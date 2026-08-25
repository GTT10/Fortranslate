# Decision 0089: plan reactive EB patches at outflow boundaries

## Context

Configured single patches could meet an outflow physical side, but the
temperature tagger skipped every root boundary cell, plan validation required
strictly internal tag and patch bounds, and patch-set hydro rejected any child
that reached the domain edge. Consequently a dynamic regrid could neither
create nor retain the physical-side topology already supported by the fine
hydro closure.

## Decision

Evaluate temperature jumps on every active root cell. At a boundary, ignore
cardinal neighbor locations outside the domain and compare only existing active
neighbors. Continue excluding covered cells and covered neighbors. Permit tag
and patch bounds over the full root index range, clamp buffer expansion to that
range, and grow minimum-size intervals toward whichever in-domain side remains
available.

Apply the same rules to deterministic multipatch flood fill and component
growth. Preserve the existing two-cell separation/coalescing contract. Remove
the patch-set hydro's physical-side rejection and reuse the current child state
for its outflow exterior closure, while flux registers continue to skip sides
without an uncovered coarse neighbor.

## Consequences

Temperature features at or near an outflow boundary can dynamically create,
move, resize, and remove single or separated fine rectangles. Boundary tags
are explicit planner inputs rather than being silently discarded. Existing
transactional topology replacement, overlap retention, conservation,
synchronization, and rollback apply without a second representation.

The one-sided tag metric is not a periodic seam treatment and does not add
inflow, wall, or mixed boundary refinement. Deeper EB levels, EB molecular
transport, distributed ownership, and explicit physical-side checkpoint parity
remain outside this decision.
