# Decision 0190: Persist the patch-tree conservation baseline

## Context

Patch-tree restart restores numerical fields and cumulative clock metadata,
but both public drivers recompute `initial_integrals` from the restored state.
Their final conservation diagnostic therefore describes only the continuation
suffix, not the complete logical run represented by the checkpoint.

## Decision

Store one finite composite-integral value for every conserved component. Move
the base checkpoint envelope from schema 2 to 3 and the public fingerprinted
envelope from schema 5 to 6. Public drivers pass the original run baseline on
every write and restore it before continuation. Low-level callers that omit the
optional value receive a baseline computed from the tree being written.

Sparse writes require communicator-wide agreement on presence, vector size,
and every value before gathering fields. The selected I/O root serializes the
vector; restart broadcasts it before topology distribution and field scatter.

## Consequences

The final conservation diagnostic spans serial and changed-rank restart
boundaries. Invalid writes do not replace a checkpoint, and rejected reads
leave the optional baseline unallocated. Older strict-schema checkpoints are
intentionally incompatible; numerical state and fingerprint contents do not
otherwise change.
