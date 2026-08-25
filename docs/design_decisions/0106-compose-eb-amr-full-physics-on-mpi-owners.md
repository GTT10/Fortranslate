# Decision 0106: compose reactive EB AMR full physics on MPI owners

## Context

Chemistry, hydrodynamics, and molecular transport have independent owner-only
MPI transactions. A usable reactive interval must compose them symmetrically
without exposing a partially advanced hierarchy or misleading counters when a
later operator rejects.

## Decision

Wrap the existing owner operators in one `R-T-H-T-R` transaction: reaction for
half an interval, SSPRK2 transport for half an interval, hydro for the complete
interval, transport for the second half, and reaction for the final half. Run
every nested transaction on private candidate root and patch-set fields.

Publish the final hierarchy, the minimum transport limiter across both half
steps, and the sums of committed chemistry, hydro, and transport owner calls
only after the final chemistry transaction succeeds. Initialize every public
counter to zero so failure at any stage has one unambiguous result.

## Consequences

The MPI EB AMR bridge now owns a complete reactive split interval with serial
multipatch parity and exact all-operator accounting on one, two, four, and
eight ranks. Invalid hydro controls after successful reaction and transport
prefixes leave caller fields bitwise unchanged and report no committed work.

The transaction still operates through replicated complete EB fields and does
not select a timestep, advance a public clock, regrid, checkpoint, or write
distributed output. Those responsibilities remain future milestones.
