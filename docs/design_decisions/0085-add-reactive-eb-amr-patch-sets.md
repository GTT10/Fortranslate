# Decision 0085: add reactive EB AMR patch sets

## Context

The reactive EB AMR lifecycle represented at most one rectangular fine patch.
A single bounding box around disconnected temperature features refines the
space between them unnecessarily, while advancing independently constructed
rectangles without a set-wide transaction risks double-counted root cells,
overlapping redistribution neighborhoods, repeated root advances, or partial
hierarchy commits.

## Decision

Add a deterministic collection planner and a separate two-level reactive EB
patch-set type. Cluster disconnected tags with a configurable gap, expand each
cluster independently, and coalesce candidates that cannot maintain the
two-coarse-cell separation needed by the 3-by-3 EB redistribution stencil.
Reject collections that overlap or touch a physical boundary.

Make creation, movement, repartition, removal, average-down, and composite
integration set-wide operations. Average every old child into a private root
before constructing the replacement set by PCM, then copy exact aligned fine
data across every old/new patch intersection. Commit only after all geometry,
state, temperature, and EOS checks succeed.

Advance the root once per coarse interval and subcycle every child from the
same old/new root interval. Give each child its own EB flux register, apply
reflux in deterministic order, and average down the complete set. Compose
active-cell chemistry around that hydro transaction with whole-hierarchy
rollback.

Keep the existing public single-patch lifecycle and checkpoint schema unchanged
until the patch-set kernel is independently qualified.

## Consequences

Disconnected features can now be represented by multiple conservative fine
rectangles without advancing the root repeatedly. Topology and time advancement
have explicit set-wide atomicity, and the separation contract keeps independent
reflux/redistribution corrections from coupling two patch interfaces.

The public application does not yet select this representation, and its
checkpoint format cannot serialize it. Deeper levels, EB molecular transport,
and distributed patch ownership also remain outside this decision.
