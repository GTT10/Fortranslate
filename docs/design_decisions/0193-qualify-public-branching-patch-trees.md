# Decision 0193: Qualify public branching patch trees

## Context

Arbitrary-depth serial and sparse-MPI libraries exercised branching trees in
direct tests, but every public full-physics and restart case formed only one
parent-child chain. The application boundary therefore did not prove that
multiple sibling identities survive scheduled regridding, sparse ownership,
formatted checkpointing, and changed-rank restart together.

## Decision

Use two separated temperature features in the established public boundary and
restart cases. Keep one feature on the x-upper physical side and one in the
interior. Require four populated levels, physical-side contact at every level,
and at least two leaf-visible patch identities on one level.

Reuse the existing reacting full-transport fresh and split-run gates. Compare
serial and 1/2/4/8-rank fresh output, independent serial restart, and two-rank
checkpoint continuation at four and eight ranks. Keep checkpoint schema 8;
ordered branching topology is already self-describing state.

## Consequences

The public lifecycle now qualifies sibling branches across process and owner
boundaries instead of relying only on direct library tests. The case still does
not represent periodic-seam patches or multiple physical boundary types.
