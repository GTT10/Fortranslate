# Decision 0093: compose three-level EB Strang chemistry

## Context

The recursive three-level kernel advanced reactive hydrodynamics and closed a
finest interface cut by the embedded boundary, but chemistry composition was
available only for root-only, two-level single-patch, and two-level patch-set
lifecycles. Calling those two-level operations independently would repeat the
root update and would not provide one transaction across the nested hierarchy.

## Decision

Add one static three-level driver operation that owns a private candidate for
every state and temperature field. Apply a reaction half-step to active cells
on all three levels, run the existing recursive three-level hydro operation,
then apply the second reaction half-step to all levels.

After the second reaction stage, synchronize the reactive hierarchy from the
finest level to the middle and then from the middle to the root. This removes
the independently reacted values in parent cells covered by a child and
recovers EOS-consistent parent temperatures. Publish all six fields only when
every chemistry, hydro, conservation-closure, and synchronization stage
succeeds.

## Consequences

The qualified static hierarchy can combine full seven-species elementary
chemistry with recursively subcycled EB hydrodynamics, including a cut finest
interface. Composite mass, total energy, density/species closure, positive
temperature, deepest-first synchronization, and exact rollback are verified.

The operation does not introduce lifecycle ownership. A public timestep loop,
three-level CFL reduction, regridding, checkpoint/restart, output, molecular
transport, arbitrary depth, and MPI distribution remain separate work.
