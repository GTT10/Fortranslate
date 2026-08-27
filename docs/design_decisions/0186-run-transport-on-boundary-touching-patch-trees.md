# Decision 0186: Run transport on boundary-touching patch trees

## Context

The public arbitrary-depth x-upper cases qualify hydrodynamic advancement,
sparse ownership, and checkpoint/restart. Recursive molecular transport is
qualified on interior library trees, but the public boundary cases leave it
disabled and therefore do not exercise physical-side diffusive context or
register omission.

## Decision

Enable thermal conduction in the fresh and split-run boundary cases. Reuse the
existing recursive SSPRK2 transport schedule, current-fine physical exterior,
coarse-time context elsewhere, `r^2` child subcycling, diffusive registers,
reflux, and average-down. Require serial and 1/2/4/8-rank fresh parity plus
serial and cross-rank checkpoint continuation.

Do not add duplicate cases. The transport-active boundary topology replaces
the hydro-only public case while still executing the complete hydro schedule.

## Consequences

The public arbitrary-depth outflow-boundary lifecycle now exercises an active
molecular transport operator across fresh execution and restart. The focused
case isolates Fourier conduction; viscous and species-diffusive boundary-tree
applications remain later qualification work.
