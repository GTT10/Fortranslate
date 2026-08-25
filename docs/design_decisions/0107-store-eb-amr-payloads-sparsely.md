# Decision 0107: store MPI reactive EB AMR payloads sparsely

## Context

Owner-only EB operators no longer duplicate physics work, but every rank still
retains complete root and child numerical fields. That correctness bridge does
not provide the memory scaling expected from distributed AMR.

## Decision

Separate replicated geometry, topology, and ownership metadata from numerical
payloads. Store one compact entry for every root tile and child, but allocate
its conserved state and temperature only when the current rank owns it. Define
the sparse object's validity against the authoritative distribution and patch
topology, including exact shape, finiteness, and allocation ownership.

Provide an owner scatter that discards nonowner replicas and an explicit
materialization transaction that reconstructs complete temporary fields by
broadcasting each owner payload. Validate communicator-wide before publishing
the reconstructed hierarchy. Report exact locally stored scalar counts for
memory-accounting gates.

## Consequences

Persistent 2D EB numerical storage can now be genuinely rank-local while
retaining a compatibility boundary for existing replicated operators. Invalid
or missing owner payloads reject without changing materialization outputs.

Physics still materializes complete temporary fields. Direct sparse physics,
point-to-point transfers, owner migration, topology changes, distributed
checkpoint/output, and a public MPI time loop remain future work.
