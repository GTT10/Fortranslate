# Decision 0057: stream sparse physics point to point

## Context

After owner migration, adjacent halos, average-down, synchronization, and final
ghost refresh became point to point, recursive sparse hydro and transport still
broadcast each patch's interval states, complete face-flux array, and cumulative
level counter. The complete flux was replicated mainly to keep flux registers
identical even though only the patch owner consumes the register at reflux.

## Decision

Keep each patch's complete flux and coarse/fine flux register authoritative on
that patch owner. Send packed interval start/end states once to every distinct
remote child owner. Each child invocation sends only its two time-integrated
boundary-flux vectors to the parent owner.

Compute adjacent shared-face fluxes on the parent owner. Apply a correction
locally when the parent owns the affected child, otherwise send only the
correction vector to that child owner. Accumulate level-counter deltas locally
during recursion and sum the complete counter array once at the physics-stage
boundary.

Expose local counts for interval-state, boundary-flux, and correction transfers.
An independent hierarchy traversal applies the same hyperbolic `r` or parabolic
`r^2` invocation weights and must match their communicator sums exactly.

## Consequences

Sparse chemistry, hydro, transport, and full split physics no longer replicate
payloads with `MPI_Bcast`. Communication follows state and flux ownership while
existing collective acceptance reductions retain transactional failure
boundaries.

Topology-changing regrid still materializes a temporary complete tree and
remains the next distributed-memory scalability boundary.
