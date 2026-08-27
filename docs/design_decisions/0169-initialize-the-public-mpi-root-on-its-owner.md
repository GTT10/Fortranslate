# Decision 0169: Initialize the public MPI root on its owner

## Context

The public sparse-MPI arbitrary-depth EB application converted to exclusive
owner storage before any child existed, but every rank still constructed a
complete numerical root field first. That temporary field replication was the
last numerical exception in fresh application startup.

## Decision

Build the geometry-only root topology and deterministic distribution first.
Call the reactive root initializer only on the selected root-node owner. Add a
collective root-only sparse initializer that requires unallocated field inputs
on non-owners, validates the owner's complete fields, and uses `move_alloc` to
transfer them into the sparse root node. Require exactly one initializer rank
and no allocated source fields after the transfer.

## Consequences

Fresh public startup no longer constructs, copies, or broadcasts a replicated
numerical root. The one owning rank temporarily holds the root arrays and then
hands their allocation directly to sparse storage. EB geometry and tree
relations remain replicated as compact routing metadata. The compatibility
initializer from an already replicated numerical tree remains available to
verification and legacy callers but is no longer used by the public app.
