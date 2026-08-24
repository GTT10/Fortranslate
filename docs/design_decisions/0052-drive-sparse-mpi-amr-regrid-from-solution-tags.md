# Decision 0052: drive sparse MPI AMR regrid from solution tags

## Context

The sparse topology transaction accepts an explicit plan, but a running AMR
application must derive that plan from the current solution. The serial planner
already performs parent-local gradient tagging, clustering, and recursive
prolongation through the configured maximum level.

## Decision

Expose one tag-driven sparse regrid entry point. Materialize the current sparse
owner state, apply the qualified serial tagged rebuild, and collectively agree
on the tagged-cell count in addition to the topology-change and overlap counts.
Commit through the same helper used by explicit-plan sparse regrid so the
solution and deterministic owner distribution remain one transaction.

An unchanged tag-derived plan updates only evaluation bookkeeping. Any invalid
tagging configuration, inconsistent count, or commit failure restores the old
sparse solution and returns the old distribution with zero reported counts.

## Consequences

Solution-driven disconnected features can now build arbitrary-depth sparse MPI
AMR topologies without an externally prepared plan. The result returns to
globally single-copy owner storage and matches serial tagging, clustering,
rebuild, and conservation behavior.

Planning still uses a temporary full reactive replica. Owner-local tagging,
distributed plan assembly, point-to-point overlap transfer, and removal of the
temporary replica remain separate work.
