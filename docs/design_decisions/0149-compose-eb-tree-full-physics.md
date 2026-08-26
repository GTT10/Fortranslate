# Decision 0149: Compose EB patch-tree full physics

## Context

The numerical EB tree independently owns arbitrary-depth active-cell
chemistry, recursive hydrodynamics, and recursive SSPRK2 molecular transport.
Calling those public transactions separately would commit intermediate fields
and counters, so a late hydro, transport, or chemistry rejection could not
roll back the complete split interval.

## Decision

Compose one private tree candidate in the established second-order order:

```text
chemistry(dt/2)
SSPRK2 molecular transport(dt/2)
recursive hydrodynamics(dt)
SSPRK2 molecular transport(dt/2)
chemistry(dt/2)
```

Reuse the qualified all-node chemistry traversal, recursive transport
transaction, and recursive hydro transaction without duplicating their
numerical kernels. Accumulate chemistry, transport-Euler, and hydro node counts
privately. Publish the complete tree, the minimum transport limiter theta, and
all three count vectors only after final deepest-first synchronization and
validation.

## Consequences

Arbitrary-depth and branching EB trees now execute one atomic `R-T-H-T-R`
full-physics interval. A runtime three-level chain retains field and
temperature parity with the established fixed-depth implementation. A separate
four-level branching gate verifies the exact chemistry, transport, and hydro
schedules plus composite conservation and thermodynamic validity. A hydro
failure after valid reaction and transport prefixes restores the accepted tree
and zeroes every public count.

This decision does not add public time/step ownership, dynamic tagging,
checkpoint I/O, or MPI ownership for the numerical tree.
