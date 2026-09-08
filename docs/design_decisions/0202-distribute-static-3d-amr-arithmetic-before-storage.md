# Decision 0202: distribute static 3D AMR arithmetic before storage

## Context

Milestone `0.208.0` made the static two-level 3D AMR hierarchy restartable,
but every hydrodynamic operation was still serial. Moving immediately to
rank-local hierarchy storage would combine decomposition, halo exchange,
coarse/fine synchronization, I/O, and numerical-order changes in one step.
That would make a parity failure difficult to localize.

## Decision

Assign contiguous x-plane slabs independently on the coarse and fine levels.
Each rank evaluates CFL rates, all directional faces anchored on its planes,
and SSPRK2 candidates for its cells only. Reconstruct complete face and state
arrays with fixed-order collectives, then run the already qualified serial
reflux and average-down sequence redundantly on each replicated hierarchy.

Before distributed arithmetic, require exact consensus over the patch,
floating controls, Riemann solver, and complete NASA7 species names and
coefficients. Before hierarchy broadcast, establish one common root and one
common fixed-size layout descriptor so differing array counts are rejected
without entering an unsafe payload collective. Receive all payloads into
private candidates and publish only after collective validation.

Keep formatted checkpoint and CSV I/O on rank zero. Checkpoints remain
owner-map-free: rank zero may read a checkpoint written by another valid rank
count and transactionally broadcast it to the new slab distribution.

## Consequences

The public static hydro case can prove serial and 1/2/4/8-rank byte identity,
including a two-rank checkpoint resumed at four and eight ranks. Failure and
metadata-consensus behavior is isolated from later storage work. The cost is
that every rank still stores both complete levels and collective reconstruction
traffic scales with the full fields. This milestone is therefore a
distributed-arithmetic reproducibility boundary, not a scalable-memory or
production-I/O claim.
