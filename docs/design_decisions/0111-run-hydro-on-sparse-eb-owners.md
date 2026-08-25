# Decision 0111: run hydro on sparse EB owners

## Context

The sparse full-physics transaction still materializes every fine child before
hydro, although child subcycling and reflux are exclusive-owner operations.
Root EB StateRedist, by contrast, uses overlapping level neighborhoods and
remains a level-wide physics operation.

## Decision

Assemble only the distributed root tiles into a level-wide temporary, advance
that root on its physics owner, and synchronize its updated state and fluxes.
Advance each fine child directly from its sparse owner allocation, including
coarse-time exterior construction, ratio subcycling, flux-register
accumulation, and reflux. Do not broadcast updated child state.

After all children succeed, copy corrected root rows back to their sparse
owners and invoke direct sparse reactive average-down. Retain the caller's
sparse object until the complete operation is accepted collectively.

## Consequences

Fine numerical state stays globally single-copy during hydro, with serial
field parity and exact rollback after a late child failure. Root arrays and
fluxes remain temporarily replicated, reflecting the current undecomposed
root StateRedist algorithm.

Root halo/redistribution decomposition, direct sparse transport, targeted
coarse/fine traffic, time-loop control, regridding, checkpointing, and output
remain future work.
