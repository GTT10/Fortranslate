# Decision 0036: react every tree patch before deepest-to-root synchronization

## Context

Reactive patch-tree hydro advances levels on different subcycle schedules, but
cell-local chemistry is not constrained by the mesh spacing and existing AMR
paths apply the same physical reaction interval to all represented cells.
Reacting only leaves would leave uncovered parent cells unchanged; reacting
only parents would discard fine-scale thermodynamic states.

## Decision

Apply each chemistry half-step independently to every patch at every level for
the same physical interval. After all patches react, extract their interiors,
average down every relation deepest-to-root, copy the synchronized fields back,
recover temperatures, and refill ghosts. Compose this operator symmetrically
around recursive hydro as `R(dt/2)-H(dt)-R(dt/2)`.

Deep-copy the complete solution before the first reaction half-step. Restore
that copy if chemistry, synchronization, hydro, temperature recovery, or ghost
fill fails at any point.

## Consequences

Uncovered parent cells and all fine cells receive chemistry while covered
parent representations remain consistent with their children between split
operators. Time and level-advance counters remain owned by hydro. This slice
uses the existing elementary chemistry integrator; transport, dynamic tree
rebuilds, same-level exchange, and distributed patch ownership remain future
work.
