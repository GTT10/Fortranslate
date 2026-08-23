# Decision 0032: rebuild a two-level patch collection from synchronized tags

## Context

The qualified multipatch engine could advance, synchronize, and conservatively
regrid a supplied set of separated patches, but the public reactive AMR
application still constructed one patch per level. Dynamic application use
requires tag clustering, empty-set handling, old/new overlap transfer, and
composite output without weakening the existing arbitrary-depth path.

## Decision

Add an explicit `amr_multipatch_enabled` configuration switch for a two-level
patch-set application. Tag the synchronized root, cluster disconnected tags
using `amr_maximum_patch_gap_cells`, expand and coalesce candidates with the
existing buffer and minimum-width rules, and constrain periodic patches away
from the seam. Outflow patches may contact either physical boundary.

Before replacing a changed collection, average every old patch down to the
root. Conservatively prolong every new patch, then copy all equal-resolution
old/new fine intersections exactly. Reevaluate the plan at the configured
coarse-step interval and permit the empty collection to remove all refinement.
Count evaluations, accepted changes, and transferred overlap cells. Emit
ordered composite output by alternating uncovered parent intervals with fine
patch rows, so each physical cell interval appears exactly once.

## Consequences

The public executable can create, move, repartition, and remove multiple
reactive patches while reusing the fixed two-level chemistry, transport, and
hydro advance. The separate arbitrary-depth path remains unchanged and owns
one patch per level. Independent adjacent boxes, same-level ghost exchange,
periodic-seam splitting, load balancing, and MPI patch ownership still require
explicit data ownership and communication designs.
