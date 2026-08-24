# Decision 0061: prescribe zero-net-mass wall species conversion fluxes

## Context

Reactive two-dimensional walls were impermeable to species. Catalytic-wall
models need a boundary transport contract before surface mechanisms, coverage
states, or reaction-rate integrators can be added. Treating wall conversion as
a total mass source would also violate the existing impermeable inviscid wall
contract.

## Decision

Give every physical face an `impermeable` or `prescribed` species mode. A
prescribed vector is expressed in kg/(m2 s), positive from the wall into the
gas, and must be finite, match the active mechanism width, and sum to zero.
Configuration rejects prescribed fluxes on periodic/inflow/outflow faces or
when molecular species transport is disabled.

Convert wall-to-gas input to coordinate-oriented flux by changing the sign at
upper faces. Keep density and momentum wall fluxes unchanged. Add
`sum(h_k J_k)` to total-energy flux using the existing face species enthalpies.
The existing face positivity limiter scales every species flux and this
species-enthalpy term by one factor.

## Consequences

The wall can convert species conservatively without adding total mass, and a
future surface-chemistry model can calculate the same boundary vector without
changing the transport update. The present interface does not compute surface
reaction rates, surface coverages, porous injection, Soret transport, or
embedded-boundary geometry.

Verification covers lower/upper orientation, zero-sum closure, enthalpy flux,
invalid-vector and disabled-transport rejection, namelist application parsing,
and a transient inventory/conservation regression.
