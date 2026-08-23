# Decision 0039: build patch-tree plans by tagging each prospective parent

## Context

The runtime patch-tree rebuild accepts arbitrary branching plans and preserves
same-resolution overlap, but callers still have to construct every
parent-local interval. A dynamic AMR engine needs to derive those intervals
from the evolving solution without making the result depend on the current
flattened patch indices.

## Decision

Average the current tree deepest-to-root and use the synchronized root as the
canonical planning state. Starting at the root, tag each prospective parent
independently with the configured normalized-gradient and absolute thresholds.
Crop the tag domain to retain parent cells required by coarse/fine ghost
interpolation, then cluster disconnected tags with the configured gap, buffer,
and minimum patch width.

Append parent-local collections in parent order and record the flattened
parent index on every child plan. Conservatively prolong the partial tree from
the root and repeat on its new children. Stop a branch when its parent has no
tags, and stop the complete recursion when no children remain or
`amr_max_levels` is reached. Pass the finished plan to the existing
transactional rebuild and physical-coordinate overlap-transfer path.

## Consequences

Planning is deterministic for a synchronized root and independent of the old
patch partition. Different parents may terminate at different depths, while
disconnected features retain separate branches. The qualified implementation
uses one configured refinement ratio for every new relation and strictly
interior children. Same-level adjacent-patch exchange, physical-boundary
children, load balancing, and MPI patch ownership remain separate.
