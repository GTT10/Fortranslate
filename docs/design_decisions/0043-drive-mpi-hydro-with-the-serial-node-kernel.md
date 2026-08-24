# Decision 0043: drive MPI hydro with the serial node kernel

## Context

Hydro couples a patch to parent time interpolation, sibling halos, shared
fine/fine fluxes, coarse/fine registers, reflux, and recursive subcycling. The
`0.49.0` ownership bridge and `0.50.0` owner chemistry establish collective
ordering and rollback, but copying the complete serial recursion into the MPI
module would create two independently evolving finite-volume algorithms.

## Decision

Extract the existing serial one-patch PCM/PLM/PPM update as a shared AMR
kernel. It returns the interval-start and interval-end states, the complete
face-flux field, and the two time-integrated boundary fluxes while incrementing
the level counter. Serial recursion calls this kernel directly.

MPI recursion walks the same deterministic patch and substep order on every
rank. Only the designated owner calls the shared kernel. After collective
acceptance, broadcast the owner's start state, face fluxes, accepted patch,
and counter. Each replica then uses the existing serial ghost-fill,
sibling-exchange, shared-flux, flux-register, reflux, average-down, and
temperature-recovery operations. Any rejected patch or collective operation
restores the owner-synchronized pre-step backup on every rank.

## Consequences

Serial and MPI hydro use one finite-volume implementation, and owner-only call
accounting is independent of communicator size. Recursive PCM and adjacent
PPM cases match serial state and conservation while exercising cross-owner
faces and deep rollback.

The hierarchy remains replicated, and non-owner ranks still apply deterministic
hierarchy and flux-register operations after owner broadcasts. This decision
does not claim sparse rank-local storage, stage-synchronous distributed Runge--
Kutta halos, owner-only molecular transport, regrid migration, or scalable
point-to-point communication.
