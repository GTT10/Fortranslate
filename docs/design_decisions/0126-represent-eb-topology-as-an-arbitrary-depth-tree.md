# Decision 0126: represent EB topology as an arbitrary-depth tree

## Context

The established 2D EB lifecycle has separate two-level multipatch and strictly
nested three-level interfaces. Those layouts qualify the numerical kernels but
cannot describe a runtime-selected depth, multiple refined parents on a deeper
level, or a complete dynamic topology replacement without adding another
fixed-depth structure.

## Decision

Introduce a geometry-only EB patch-tree topology. Store the root geometry and
a runtime-sized sequence of refinement relations. Each relation records one
refinement ratio, flattened ordered child nodes, explicit parent indices, and
parent-to-child offsets. Every node owns only its EB geometry and AMR patch box.

Construct every patch against its actual parent geometry. Require children of
the same parent to retain the established two-coarse-cell separation. Provide a
whole-tree rebuild operation that stages and validates a candidate before
commit, reports identical plans as no-ops, and leaves the accepted topology
unchanged on invalid parent links or geometry.

## Consequences

Geometry metadata can now express branching at arbitrary depth without adding
new fixed-level types. Parent/local-child lookup remains deterministic and a
future owner map can use the same flattened indexing.

This decision does not add arbitrary-depth numerical field containers,
conservative prolongation/restriction/overlap migration, physics recursion,
timestep selection, checkpoint I/O, or MPI ownership. Those operations remain
subsequent work and continue to use the qualified fixed-depth paths meanwhile.
