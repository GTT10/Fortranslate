# Decision 0163: Restart EB trees without owner metadata

## Context

The serial arbitrary-depth EB checkpoint contains complete numerical topology
and fields, but sparse MPI runs retain each node only on one deterministic
owner. Persisting that owner map would bind a restart to the writer's rank count
and work-weight exponent. Materializing the tree on every rank would also break
the sparse lifecycle boundary.

## Decision

For writes, gather each non-root-owned numerical node directly to one selected
I/O root and invoke the qualified serial writer there. For reads, let only that
root parse the checkpoint. Broadcast compact topology and complete EB geometry,
recompute the distribution from the current communicator and requested work
exponent, then scatter each numerical node directly from the root to its new
owner.

Establish collective agreement on root, lifecycle metadata, depth, work
exponent, and ordered species before file I/O. Count one transfer for every
node whose owner differs from the selected root. Publish distribution, sparse
fields, metadata, and transfer counts only after all ranks validate the result.

## Consequences

Checkpoint files contain no MPI rank count or owner metadata and can restart
under a different deterministic distribution. Only compact geometry is
replicated; conserved state and temperature never become an all-rank numerical
tree. The selected root temporarily owns one complete tree at the explicit I/O
boundary. Composite output remains separate lifecycle work.
