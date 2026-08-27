# Decision 0167: Expose the sparse MPI EB tree application

## Context

Sparse arbitrary-depth EB ownership, physics, regridding, checkpoint/restart,
and output were qualified only through a verification executable. Users could
not run that lifecycle from the public reactive 2D namelists.

## Decision

Add a separate installed `pelef_mpi_reactive_eb_patch_tree_2d` executable.
Reuse the serial patch-tree inputs and add a bounded MPI work exponent. Convert
the fresh root field to sparse ownership before recursive tagging, then use
owner-local APIs for the full lifecycle and selected root zero for file I/O.

## Consequences

The arbitrary-depth MPI path is now runnable without changing the verification
harness or fixed-depth applications. Fresh startup temporarily constructs a
replicated root field, but no recursively refined child field is replicated.
Rank-neutral application restart receives a separate end-to-end gate later.
