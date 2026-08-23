# Decision 0029: distinguish physical and coarse/fine sides per AMR relation

## Context

The original hierarchy required every fine patch to be strictly interior, so
both outer fine faces always participated in a flux register and every ghost
value came from the parent. Removing only that geometry restriction would make
the physical side index a nonexistent parent neighbor and apply an invalid
reflux correction.

## Decision

Boundary contact is derived independently for the left and right side from the
relation geometry. An outflow physical side receives fine-level constant
extrapolation for the face ghost and all four PPM/WENO stencil ghosts. Its flux
register entry is reset but never applied. Any opposite interior side keeps the
qualified parent interpolation and reflux path.

Endpoint parent cells are conservatively prolonged with zero limited slope on
the missing side. Dynamic PPM/WENO planning may remove the strict-nesting
buffer only where a parent edge coincides with a global outflow boundary.

## Consequences

A single patch may touch either or both outflow boundaries and nested children
may retain that contact. Full-parent periodic refinement can self-wrap, but a
one-sided periodic patch remains unsupported because its physical seam is also
a coarse/fine interface requiring topology-aware cross-seam reflux.
