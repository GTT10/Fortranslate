# Decision 0164: Write one composite EB tree CSV

## Context

The arbitrary-depth 2D EB tree stores complete fields on every node, including
parent cells replaced by finer children. Writing every node verbatim would
double-count refined regions, while writing one file per node would leave
consumers to reconstruct topology. Sparse MPI also retains each numerical node
on only one owner and must not materialize output fields on every rank.

## Decision

Write one deterministic CSV containing every finest-available AMR cell exactly
once. For each node, mask the union of its direct children's coarse rectangles
and write every remaining cell with level, patch, local index, spacing,
coordinate, time, EB diagnostics, conserved fields, recovered primitive state,
temperature, and ordered species mass fractions.

For sparse MPI, establish collective agreement on root, time, and species,
gather each non-root-owned node directly to the selected writer root, and call
the serial writer only there. Broadcast the file result and publish transfer
counts only after success.

## Consequences

The result is immediately consumable as one composite table without duplicate
coarse coverage. EB-covered cells remain explicit diagnostic rows. Only the
writer root temporarily allocates a complete tree, with one entity transfer
per node owned elsewhere; the owner-local physics lifecycle is unchanged.
