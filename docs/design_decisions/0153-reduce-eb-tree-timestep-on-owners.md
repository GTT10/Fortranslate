# Decision 0153: Reduce the arbitrary-depth EB timestep on owners

## Context

The MPI arbitrary-depth tree now stores numerical fields only on node owners,
but the qualified combined hydro and transport timestep selector consumes a
complete serial tree. Materializing every field before each public time step
would restore replication on the most frequently repeated control path.

## Decision

Evaluate each active node's hyperbolic and explicit-transport stability limits
only on its owner. Convert a node-local interval to the root interval using the
cumulative ancestor refinement product, take the rank-local minimum, and use a
communicator minimum for the public result. Count evaluated active nodes
separately so ranks without work can contribute a neutral huge value while an
entirely inactive tree still rejects.

Require collective agreement on the replicated distribution, sparse layout,
CFL values, species count, and enabled transport terms before evaluation.
Publish the timestep and local evaluation count only after the global minimum
is finite, positive, non-huge, and backed by at least one active node.

## Consequences

Stable-step selection no longer allocates or broadcasts a complete numerical
tree. Its numerical result remains identical to the serial all-node scan and
works after direct owner migration, including ranks that own no nodes.

Recursive hydro, transport, and chemistry still require complete-tree entry
points and remain the next owner-local routing work.
