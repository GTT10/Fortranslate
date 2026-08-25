# Decision 0101: add sibling-patch EB AMR transport

## Context

The two-level EB AMR transport transaction already subcycles one fine patch,
fills its exterior state from a time-interpolated coarse solution, refluxes a
diffusive flux register, and synchronizes the hierarchy. Multipatch EB AMR
owns a set of nonoverlapping sibling rectangles and must advance all of them
without duplicating the coarse transport update or assigning a conserved
volume to more than one child.

## Decision

Advance the coarse transport stage once and retain its old and Euler-candidate
states for every child exterior fill. Give each sibling an independent flux
register and ratio-subcycle it over the same coarse interval. Reflux each
register into a private coarse candidate; the operations act on disjoint
interface regions because the patch-set validator excludes overlapping or
touching children. Average down the complete patch set only after every child
has completed.

When any sibling interface crosses the embedded boundary, measure one global
patch-set composite residual after all reflux and average-down operations.
Apply a single density/species-consistent correction per fluid volume over
active, unrefined coarse recipients and recover their temperatures through the
EOS. This avoids applying independent conservation closures to a shared
coarse hierarchy. Compose two complete synchronized Euler transactions as
SSPRK2 and use the result in both transport half-steps of the multipatch
reactive Strang driver.

## Consequences

The public dynamic two-level multipatch lifecycle can now run the qualified
mixture molecular-transport subset with per-child subcycling, conservative
diffusive synchronization, and whole-hierarchy rollback. Its coarse timestep
includes every child's ratio-scaled parabolic limit.

Transport checkpoint/restart, same-level diffusive ghost exchange between
touching siblings, coarse-to-fine spatial slopes, thermal or catalytic
embedded walls, distributed ownership, and parallel flux registers remain
outside this milestone.
