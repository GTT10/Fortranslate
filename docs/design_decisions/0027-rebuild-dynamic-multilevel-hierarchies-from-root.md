# Decision 0027: rebuild changed multilevel hierarchies from a synchronized root

## Status

Accepted for PeleF `0.32.0`.

## Context

The recursive reactive engine can advance any strictly nested hierarchy, but
the runnable application needs to create and change that hierarchy from state
tags. Preserving every old fine value across simultaneous changes at several
depths requires a general geometric overlap mapper that does not yet exist.

## Decision

Retain the established overlap-preserving two-level path when
`amr_max_levels = 2`. For larger values, recursively tag and prolong one child
per level until tags end or the configured limit is reached. Suppress tags on
interior patch-edge cells to preserve a coarse/fine ghost neighborhood.

At regrid time, average the old hierarchy down from deepest level to root. Plan
the replacement hierarchy from that synchronized root. If its topology is
unchanged, retain the complete old hierarchy. If it changed, conservatively
prolong every child from the synchronized parent. Emit composite output with a
recursive left-parent/child/right-parent traversal.

## Consequences

- hierarchy creation, movement, depth growth, and depth reduction conserve the
  composite state;
- the old hierarchy remains untouched when its plan is unchanged;
- changed hierarchies preserve parent cell averages but can lose fine-scale
  detail even where old and new patches overlap;
- the application now runs arbitrary configured depth without changing the
  default two-level behavior;
- multiple patches, general overlap transfer, physical-boundary refinement,
  MPI patch ownership, and load balancing remain future work.
