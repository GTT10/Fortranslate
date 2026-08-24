# Decision 0041: stage MPI AMR through owner-authoritative replicas

## Context

The serial patch-tree engine owns recursive advancement, regridding, overlap
transfer, adjacent-sibling ghost fill, and fine/fine flux reconciliation. Its
state layout still allocates every patch together. Replacing that layout and
the recursive call graph in one step would mix ownership correctness, message
ordering, physics scheduling, rollback, and migration failures.

## Decision

Introduce MPI ownership outside the serial AMR modules. Require every rank to
hold an identical valid hierarchy, verified collectively from the complete
integer topology and physical root extent. Visit patches in deterministic
level/flattened order and give each patch to the rank with the least accumulated
cell work, breaking ties by rank.

During this bridge phase, retain a field replica on every rank but treat only
the assigned owner as authoritative. Synchronize complete patch interiors with
owner-rooted broadcasts. For every adjacent sibling face, broadcast the left
owner's last cells and the right owner's first cells in common hierarchy order,
producing explicit opposite-side halo layers. Reject inconsistent hierarchy
metadata collectively before any patch communication.

## Consequences

Patch identity, load accounting, owner authority, rank-crossing adjacency, and
collective ordering are independently testable for 1, 2, 4, and 8 ranks. The
serial AMR implementation remains usable without MPI, and later recursive
integration has a stable communication contract.

Replicated fields do not reduce memory use and broadcasts are not the final
scalable exchange. Owner-only hydro/chemistry/transport, stage-synchronous
shared-flux communication, sparse local storage, and patch migration after
regrid remain required before claiming fully distributed AMR execution.
