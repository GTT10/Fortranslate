# Decision 0048: advance recursive hydro on sparse MPI AMR patches

## Context

The `0.55.0` sparse chemistry path established owner-local physics,
deepest-to-root synchronization, and sparse rollback. Hydro additionally needs
time-interpolated parent states, level-ratio recursion, coarse/fine flux
registers, adjacent fine/fine flux ownership, reflux, and average-down.

## Decision

Run the existing one-patch finite-volume kernel only on the sparse owner.
Broadcast its interval-start state, accepted interval-end state, and effective
face flux for deterministic recursion. Keep compact flux-register arrays
replicated and accumulate identical coarse and fine boundary integrals on
every rank. Fill child and sibling ghosts only where their payload is local.

Reconcile adjacent time-integrated face fluxes on every rank but apply the
conservative cell correction only on each target owner. After child subcycles,
stream child interiors to the parent owner, which performs reflux,
average-down, and temperature recovery. Preserve a full local sparse backup
and require collective acceptance at every mutation boundary.

## Consequences

Recursive hydro no longer materializes a complete reactive tree. Mixed-ratio
subcycling, adjacent PPM children, exact counters, and rollback retain serial
behavior while persistent payloads remain owner-local.

Interval data and child interiors are still sent through collectives, and
flux-register metadata remains replicated. Molecular transport, combined full
physics, point-to-point schedules, and topology-changing distributed regrid
remain separate work.
