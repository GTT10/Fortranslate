# Decision 0151: Assign MPI owners before sparse tree storage

## Context

The arbitrary-depth EB numerical tree now owns complete serial physics and a
public clock, while the established MPI EB path is limited to fixed-depth
patch sets. Moving recursive physics directly to sparse storage would combine
owner selection, memory layout, migration, context routing, and numerical
correctness in one change.

## Decision

Introduce a topology-matched distribution first. Assign every runtime node to
one deterministic rank by greedy accumulated work. Count allocated cells as
the base cost and optionally weight deeper levels by the first or second power
of their cumulative subcycle product.

Keep a complete numerical tree on every rank for this milestone. Treat the
assigned owner as authoritative for each node, broadcast state followed by
temperature into a private replicated candidate, and commit only after
collective input and final validation. Publish exact per-rank cell, node, work,
and local owner-publication counts.

## Consequences

Arbitrary-depth and branching topologies now have a qualified MPI ownership
contract at one, two, four, and eight ranks. Control disagreement or one
rank's invalid candidate rejects collectively without exposing partial
publication. The owner map and accounting can be reused by sparse allocation
and direct migration.

Nonowners still allocate complete fields, and broadcasts still replicate every
node. Sparse owner storage, point-to-point migration, distributed timestep
selection, and owner-local recursive physics remain follow-on work.
