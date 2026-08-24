# Decision 0049: advance molecular transport on sparse MPI AMR patches

## Context

The sparse chemistry and hydro paths established owner-local physics,
distributed AMR synchronization, recursive subcycling, and rollback. Molecular
transport has the same recursive topology as hydro but uses an SSPRK2
diffusive kernel and `r²` child subcycling for each refinement relation.

## Decision

Run each transport patch interval only on its sparse owner. Stream the
interval-start state, accepted interval-end state, and effective diffusive
flux. Execute `r²` child intervals in serial order with time-interpolated
parent ghosts. Keep compact diffusive registers replicated, reconcile adjacent
time-integrated diffusive faces, and apply state corrections only on owners.

After each child schedule, stream child interiors to the parent owner for
diffusive reflux, average-down, and temperature recovery. Reuse sparse final
ghost refresh, collective acceptance, local call accounting, and exact sparse
rollback.

## Consequences

All three component physics operators now run directly on rank-local AMR
payloads with qualified recursive synchronization and failure behavior.

The combined full-physics transaction still uses replicas. Interval data,
fluxes, and compact registers retain the correctness-first collective schedule;
point-to-point communication and topology-changing distributed regrid remain
separate work.
