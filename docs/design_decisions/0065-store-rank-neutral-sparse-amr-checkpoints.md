# Decision 0065: store rank-neutral sparse AMR checkpoints

## Context

The public sparse MPI AMR driver can evolve and regrid owner-local patch
payloads, but a process interruption previously required restarting from the
initial condition. Persisting the current owner map would couple a checkpoint
to the communicator size that wrote it and would prevent flexible queue or
resource changes at restart.

## Decision

Define a versioned formatted checkpoint for the complete reactive patch tree.
Store mechanism species names and conserved width, root geometry, refinement
plans, interior conserved states and temperatures, physical time, coarse-step
count, regrid/overlap accounting, and per-level hydro and transport advances.
Do not store MPI ranks or patch owners.

At a scheduled cadence, collectively materialize the authoritative sparse
tree and let rank zero replace the configured checkpoint file. A clean
stop-after-write option supports split-run regression and batch-job chaining.
On restart, validate the schema, mechanism layout, root geometry, and supported
hierarchy depth, reconstruct the replicated tree, recover ghosts and
temperatures, then run the deterministic ownership scheduler for the current
communicator before scattering owner-local payloads.

## Consequences

The numerical state and AMR accounting survive restart while ownership can
change from two ranks to four or eight. Normal timesteps and regrids remain
field-sparse; only checkpoint boundaries materialize a complete tree.

The schema is intentionally implementation-specific rather than an AMReX
plotfile format. Rank-zero formatted output is portable and easy to validate,
but is not scalable parallel I/O and does not provide an atomic multi-file or
redundant failure-recovery protocol.
