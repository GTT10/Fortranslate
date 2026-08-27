# Decision 0156: Route EB tree hydro on owners

## Context

An arbitrary-depth hydro node needs more than its own state. Each refined edge
also needs time-interpolated parent boundary context, coarse and fine flux
registration, bidirectional reflux updates, ordered average-down, and subtree
conservation closure. Materializing the complete tree would defeat sparse
ownership at every substep.

## Decision

Run the serial recursive schedule collectively, but execute each numerical
kernel only on the node owner. Extract and transfer compact parent start/end
exterior context once per distinct-owner edge. Return fine fluxes to the parent
owner after every child substep, keep the flux register there, and apply reflux
to the parent-authoritative state. Send current child state to the parent for
reflux, return the corrected child, and send it once more for ordered average-
down. Keep shared-owner edges entirely local.

Measure conservation with prevalidated owner-local subtree reductions and
apply the established closure only on unrefined active cells of the parent
owner. Traverse every rank through identical communication and acceptance
boundaries. Commit the sparse candidate and all counters only after the root
recursion validates.

## Consequences

Recursive hydrodynamics no longer materializes numerical nodes on nonowners.
A distinct-owner edge costs `r + 4` grouped direct payloads per parent
invocation; shared ownership costs none. The result follows the serial subcycle,
reflux, restriction, and conservation order within qualified floating-point
roundoff. Molecular transport still needs its separate SSPRK2 flux and limiter
routing milestone.
