# Decision 0104: run reactive EB AMR hydro on MPI owners

## Context

Reactive chemistry now runs on EB AMR owners, but hydrodynamics still runs
only through the serial multipatch transaction. A root y-tile is not yet a
safe independent hydro entity: weighted StateRedist forms overlapping level
neighborhoods and its second-order reconstruction may inspect cells beyond a
tile edge. Treating those edges as physical boundaries would change both the
finite-volume flux and conservative redistribution.

## Decision

Keep root storage synchronization tiled, but select the owner of the first
root tile as the exclusive root-level physics owner. That rank advances the
complete root EB level once and publishes its state, recovered temperature,
and face-centroid fluxes. Treat each fine sibling as an independent owner
physics entity. Its owner builds the coarse flux contribution, performs all
refinement-ratio substeps with time-interpolated root exterior data,
accumulates the fine flux contribution, and refluxes the current root
candidate. Child results are committed in deterministic patch order. The root
physics owner performs the final set-wide average-down.

Validate all numerical and string controls collectively before execution.
Retain the caller's root and child fields until the root, every child, reflux,
and final average-down have all been accepted communicator-wide. Report only
committed level advances.

## Consequences

The root finite-volume/StateRedist transaction and every fine subcycle now
execute exactly once on exclusive MPI owners. The operation preserves the
serial multipatch order and reaches serial field parity on one, two, four, and
eight ranks. A late fine-owner failure discards the already computed root and
earlier-child candidates on every rank.

The root level is intentionally one physics entity, so this milestone removes
replicated hydro computation without yet decomposing weighted StateRedist.
Sparse EB storage, point-to-point root halos, decomposed StateRedist,
owner-only molecular transport, topology migration, and a public distributed
time loop remain future work.
