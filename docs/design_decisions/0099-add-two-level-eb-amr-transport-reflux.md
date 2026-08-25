# Decision 0099: add two-level EB AMR transport reflux

## Context

The single-level EB transport operator already owns mixture transport,
face-centroid interpolation, fluid-inventory limiting, StateRedist, and EOS
recovery. A two-level hierarchy also needs fine exterior data and conservative
coarse/fine diffusive synchronization. Duplicating the transport physics in
the AMR layer or refluxing only hydrodynamic fluxes would create inconsistent
operators and lose composite conservation.

## Decision

Expose the single-level transport face-centroid fluxes together with its
conservative right-hand side. Represent spatially varying transport exterior
data independently of physical boundary kinds, and recover it from the
existing EB coarse-to-fine conserved exterior state.

Advance one coarse transport Euler interval and `r` fine Euler substeps using
time-interpolated coarse exterior states. Accumulate the coarse and fine
diffusive fluxes in the established EB flux register, then apply reactive
reflux and average-down. Compose two complete synchronized hierarchy Euler
transactions as SSPRK2. Select the root interval from both the coarse limit and
`r` times the fine limit, and insert transport half steps into the existing
reactive Strang driver.

## Consequences

The public single-patch two-level EB AMR lifecycle now preserves composite
conserved quantities across molecular-transport interfaces and retains
transactional failure behavior. Physical outflow sides use the established
zero-gradient fine exterior, while interior patch sides use time-interpolated
coarse data.

The milestone does not add coarse-to-fine spatial slopes, thermal or catalytic
embedded walls, sibling-patch transport, three-level transport, distributed
ownership, transport checkpoint/restart, or parallel flux registers.
