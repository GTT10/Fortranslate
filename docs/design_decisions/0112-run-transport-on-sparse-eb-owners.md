# Decision 0112: run transport on sparse EB owners

## Context

The sparse full-physics transaction still materializes every fine child before
molecular transport, although each child's SSPRK2 stages, ratio subcycles, and
diffusive reflux are exclusive-owner operations. Root transport and EB
StateRedist still require overlapping level-wide neighborhoods.

## Decision

For each SSPRK2 Euler stage, assemble only the distributed root tiles into a
level-wide temporary. Advance root transport and StateRedist on the root physics
owner and synchronize the updated root state, temperature, and fluxes.

Advance each fine child directly from its sparse owner allocation. Keep
coarse-time exterior construction, ratio subcycling, diffusive flux-register
accumulation, reflux, and temperature recovery on that owner; do not broadcast
the child payload. Synchronize the corrected root between children, then return
its rows to sparse tiles and average down directly.

For an EB-cut coarse/fine interface, compute the composite integral from local
root tiles and child payloads with a communicator reduction. Apply the closure
only to uncovered, unrefined root cells. Retain the caller's sparse object and
publish the limiter and Euler count only after both stages and the final blend
are accepted collectively.

## Consequences

Fine transport state stays globally single-copy throughout SSPRK2, with serial
state, temperature, and limiter parity and exact rollback after a late child
failure. Root arrays and fluxes remain temporarily replicated, and composite
closure introduces a reduction whose ordering may differ from serial summation.

Root transport/StateRedist decomposition, targeted root and coarse/fine
traffic, use of this direct operator in the outer sparse full-physics
transaction, public time-loop control, regridding, checkpointing, and output
remain future work.
