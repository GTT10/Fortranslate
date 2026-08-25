# Decision 0109: average down sparse EB owners directly

## Context

Chemistry executes on sparse root tiles and children, but its synchronization
materializes the complete numerical hierarchy on every rank. Reactive
average-down needs only a child's fine state, its geometry, and the intersecting
coarse owner cells.

## Decision

Let each child owner compute the volume-fraction-weighted conserved state over
the child's coarse footprint in the same component and cell order as the serial
kernel. Broadcast that compact restriction buffer, preserve covered coarse
cells, and let intersecting root tile owners recover temperature from their
local coarse guesses before applying the restricted state.

Keep a deep sparse backup and candidate until every child restriction and EOS
recovery is accepted communicator-wide. Chemistry invokes this sparse
average-down directly and never crosses the complete materialization boundary.

## Consequences

Persistent and chemistry-temporary state remains owner-local except for one
coarse-footprint conserved buffer per child. The result is bitwise-identical to
serial reactive average-down, and a late EOS failure restores every sparse
allocation exactly.

The compact buffer is currently broadcast to the communicator. Targeted
point-to-point child/root traffic, sparse hydro and transport, distributed time
advancement, topology migration, checkpointing, and output remain future work.
