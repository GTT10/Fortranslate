# Decision 0168: Qualify public sparse-MPI cross-rank restart

## Context

The installed sparse-MPI arbitrary-depth EB application exposed checkpoint
and rank-neutral restart entrypoints, but its public gate covered only fresh
runs at one, two, four, and eight ranks. The application process boundary had
not demonstrated redistribution under a changed communicator and ownership
policy.

## Decision

Use the established four-level application restart case as an uninterrupted
one-rank reference, a two-rank checkpoint-stop run, and independent four- and
eight-rank restart runs. Write the checkpoint with uniform node weighting and
restart with depth-squared weighting. Extend the identity-keyed checker to
accept multiple restarted outputs and compare each complete numeric field set
with the reference.

## Consequences

The public application now qualifies root-only checkpoint read, owner-map
recomputation, direct sparse field scatter, resumed regridding and physics,
and deterministic composite output across both rank-count and ownership-policy
changes. Fresh startup still has one temporary replicated root field; removing
that allocation remains separate work.
