# Decision 0184: Qualify arbitrary-depth outflow-boundary children

## Context

The shared 2D EB tag planner can clip a child rectangle to the parent domain,
and the recursive hydro path can already fill a physical child side from its
current fine boundary state. The public arbitrary-depth cases exercised only
strictly interior children, so neither the serial lifecycle nor sparse owner
routing had a regression contract for this topology.

## Decision

Add a four-level public hotspot case whose descendants meet the x-upper
outflow boundary. Require every populated level to reach that side, retain
coarse-time exterior interpolation elsewhere, omit physical-side register and
reflux work, and compare the sparse-MPI composite output at one, two, four,
and eight ranks.

Reuse the existing domain-inclusive planner, exterior builder, flux register,
and direct sparse ownership routes. Do not add a separate boundary topology or
silently reinterpret non-outflow boundary conditions.

## Consequences

Outflow-boundary children become a qualified public arbitrary-depth topology
without introducing a second execution path. Non-outflow refined boundaries,
periodic-seam children, and checkpoint continuation of a boundary-touching
tree remain separate milestones.
