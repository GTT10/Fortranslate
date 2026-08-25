# Decision 0105: run reactive EB AMR transport on MPI owners

## Context

Reactive EB AMR chemistry and hydrodynamics execute on exclusive MPI owners,
but the molecular-transport transaction remains serial. Transport has the same
level-wide weighted-StateRedist constraint as hydro and additionally requires
two synchronized SSPRK2 Euler stages, diffusive flux-register reflux, and one
global conservation closure when sibling interfaces cross the embedded
boundary.

## Decision

Treat the complete root level as one transport physics entity owned by the
first root-tile owner. In each Euler stage that owner computes mixture
viscosity, thermal conduction, species diffusion, barodiffusion, the EB-aware
right-hand side, StateRedist, and root face fluxes exactly once. Each child
owner builds time-interpolated root exterior data, performs every ratio
substep, accumulates an independent coarse/fine diffusive flux register, and
refluxes the current root candidate. Publish child results in deterministic
patch order.

After all children succeed, let the root owner average down the patch set and
apply the existing set-wide EB-cut conservation closure. Execute two complete
synchronized Euler transactions, blend the original owner-authoritative state
with the second Euler result on the relevant owners, recover temperature, and
average down once more. Validate transport records, boundary payloads,
switches, interval, and StateRedist controls collectively before execution.
Do not publish caller fields or advance counts until the entire SSPRK2
transaction succeeds.

## Consequences

Qualified mixture molecular transport now executes once per root or child
interval on exclusive MPI owners and retains serial multipatch field and
limiter parity at one, two, four, and eight ranks. Stale nonowner replicas do
not affect the SSPRK2 blend, and a late child failure discards all earlier
root, child, reflux, and stage candidates communicator-wide.

The correctness bridge still broadcasts complete fields. Sparse EB storage,
point-to-point root/child traffic, decomposed weighted StateRedist, distributed
flux registers, dynamic topology migration, parallel checkpoint I/O, and a
public MPI EB AMR time loop remain future work.
