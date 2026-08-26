# Decision 0127: migrate reactive state through the EB patch tree

## Context

The geometry-only patch tree can describe arbitrary depth and branching, but
the accepted reactive fields still live in fixed two-level multipatch or
strictly nested three-level structures. Dynamic tree replacement therefore
needs a numerical representation and a conservation rule before runtime
physics or distributed ownership can use the new topology.

## Decision

Mirror every topology level with an ordered numerical level. Store conserved
state and recovered temperature in each numerical node while retaining all EB
geometry and patch metadata in the topology.

Initialize children from their actual parents with the established reactive EB
PCM operation. Synchronize a hierarchy by restricting every relation in
deepest-first order. Define the arbitrary-depth composite integral as the root
integral plus each child integral minus the parent region that child replaces.

Treat dynamic replacement as one transaction. Collapse a private copy of the
old tree, initialize the candidate from the collapsed root, and process each
new level in parent-first order. Retain old same-resolution cells only when
their physical grids align and their local cell and surrounding-face EB
metrics agree. Recover active temperatures with the NASA7 EOS, synchronize the
candidate, compare every conserved integral, and commit only after all checks
succeed. Identical plans are exact no-ops.

## Consequences

Arbitrary-depth branching trees can now hold, synchronize, and conservatively
migrate reactive numerical state without adding another fixed-depth API.
Failures preserve the accepted topology and fields, and a later owner map can
reuse the same level and flattened-patch indexing.

This decision does not recurse hydrodynamics, chemistry, or molecular
transport through the tree. It also does not select timesteps, tag new patches,
write checkpoints, or distribute nodes across MPI owners. Those remain
subsequent work.
