# Decision 0117: compute the sparse EB timestep on owners

## Context

The sparse MPI reactive EB hierarchy can advance every qualified physics
operator without replicating fine payloads, but it has no public stable-step
selection. A caller would otherwise need to materialize the hierarchy or
supply a manually chosen interval before entering the sparse transaction.

## Decision

Gather root tiles once to the root physics owner. Evaluate the root EB
hyperbolic CFL limit and, when enabled, the molecular-transport stability limit
only there. Evaluate the same limits for each fine child only on that child's
owner. Multiply a child's stable fine step by its refinement ratio so every
candidate represents a coarse-level interval, then select one communicator
minimum.

Require communicator consensus for CFL controls and physics flags. Reuse the
qualified transport database and boundary-control consensus before numerical
evaluation. Count packed root-tile sends, but publish dt and the transfer count
only after every owner accepts its local fields and the global minimum is
finite and positive.

## Consequences

Timestep selection no longer requires replicated fine payloads or all-rank
root arrays. The result has the same hydro/transport stability meaning as the
serial patch-set selector and is ready to drive a future clipped multi-step
sparse transaction.

The routine does not yet own time, step counters, dynamic topology,
checkpointing, or output. Those remain separate driver responsibilities.
