# Decision 0103: run reactive EB AMR chemistry on MPI owners

## Context

The two-dimensional EB AMR MPI bridge assigns every root y-tile and fine
sibling patch to one rank, but `0.110.0` only synchronizes externally modified
payloads. The first direct distributed physics operator should avoid spatial
halo and shared-flux dependencies while proving collective transaction and
rollback behavior on the existing ownership map.

## Decision

Run active-cell chemistry on the exclusive owner of each root tile and each
fine child. Keep the serial chemistry kernel unchanged: root owners call it on
their full-width y-section with the matching EB active mask, while child
owners call it on the complete patch. All ranks enter the same entity order,
collectively accept each reactor transaction, and then broadcast the accepted
owner state and recovered temperature.

Validate interval, tolerance, species count, and reaction count collectively
before execution. Retain the caller's fields until all entities have advanced
and the replicated fine-to-root average-down has succeeded. Report only
committed owner calls; any owner or synchronization failure returns the input
fields exactly.

## Consequences

Reactive source integration is no longer replicated computation. One-, two-,
four-, and eight-rank runs produce the same patch-set result as the serial
active-cell chemistry path, and a late owner rejection rolls back earlier
candidate updates on every rank.

The field payload remains replicated and broadcasts are still the
correctness-first communication bridge. Owner-only EB hydrodynamics and
molecular transport, sparse rank-local allocation, point-to-point coarse/fine
traffic, distributed flux registers, and a public distributed time loop remain
future milestones.
