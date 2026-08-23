# Decision 0033: compose parent-owned patch sets into an arbitrary-depth tree

## Context

The arbitrary-depth reactive hierarchy stores one patch per level, while the
multipatch hierarchy stores several separated children over one root. Extending
the latter by treating every level as one global patch set loses the fact that
a fine patch must be nested inside one specific parent patch and that different
parents may stop refining at different depths.

## Decision

Represent each level relation as one child patch set per flattened parent
patch. Store prefix offsets that map `(parent patch, local child)` to one
deterministic flattened child index at the next level. Every child set uses the
parent patch's local cell count and physical bounds, so existing conservative
patch-set prolongation and restriction remain reusable. Permit empty child sets
for individual parents, but require at least one child across every represented
relation.

Walk root-to-leaf for conservative prolongation and deepest-to-root for
average-down. Compute the composite integral by starting with the root and, for
every child, subtracting its covered parent interval and adding its fine-cell
integral. Validate parent ownership, physical extents, levels, spacing, and
refinement ratios before any field operation.

## Consequences

The geometry and field layer can represent branching multipatch hierarchies of
arbitrary depth with a different refinement ratio at each relation. Parent
ownership is explicit enough for later recursive time integration and MPI
partitioning. This milestone does not yet provide flux-register recursion,
reactive advancement, dynamic tree rebuilds, adjacent-box exchange, or
distributed ownership.
