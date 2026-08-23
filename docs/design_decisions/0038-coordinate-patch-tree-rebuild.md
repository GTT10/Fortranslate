# Decision 0038: transfer changed patch trees by physical cell coordinates

## Context

A runtime branching plan can move, repartition, add, or remove patches and can
change which parent owns a later-level patch. Patch indices therefore do not
identify persistent data. Reprolonging every new patch is conservative but
unnecessarily discards resolved state where old and new grids coincide.

## Decision

Synchronize the old tree deepest-to-root, build and prolong the new tree from
that root, then inspect every old/new patch pair at each common level. When the
level spacings match, intersect their physical bounds, require cell-aligned
offsets, and copy the overlapping state and temperature cells exactly. Finish
with deepest-to-root average-down.

Treat the complete rebuild as one transaction. Preserve time and physics
advance counters, record evaluations, accepted rebuilds, and transferred cell
counts, and restore the old tree on any invalid plan or failed operation.

## Consequences

Overlap survives changes in flattened indices and parent ownership without a
global patch identity system. Changed refinement ratios conservatively fall
back to root prolongation. The caller still supplies the branching plan;
automatic per-parent tagging and clustering remain the next integration.
