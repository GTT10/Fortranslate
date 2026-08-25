# Decision 0108: run chemistry directly on sparse MPI EB owners

## Context

MPI EB state can be stored sparsely, but all physics still materializes a
complete hierarchy before execution. Chemistry is cell-local and therefore the
first operator that can remove that compatibility step without halo exchange.

## Decision

Advance every locally allocated root tile and child directly through the
existing active-cell chemistry kernel. Slice the root EB mask to the owned
y-tile and use each child's local geometry mask. Validate numerical controls
and mechanism widths collectively, accept each entity communicator-wide, and
retain a complete sparse backup until all reactors and synchronization succeed.

For this milestone, materialize only after all reactions to reuse the qualified
reactive average-down. Let the root physics owner perform that average-down,
broadcast the synchronized root, and immediately scatter root and child fields
back to a sparse candidate before publication.

## Consequences

Reactor work and persistent numerical storage are both owner-local. The result
matches serial patch-set chemistry bitwise and a late child failure restores
every sparse allocation exactly at one, two, four, and eight ranks.

Post-reaction synchronization still creates a temporary complete hierarchy.
Direct child-to-root restriction, sparse hydro and transport, time advancement,
topology migration, checkpointing, and output remain future work.
