# Decision 0098: add single-level EB molecular transport

## Context

The regular two-dimensional reactive solver already owns qualified mixture
viscosity, thermal conduction, species diffusion, and barodiffusion, while the
EB solver previously rejected all molecular transport. Reimplementing the
coefficient and face-flux physics inside EB code would create two numerical
authorities and still leave small cut cells unstable.

## Decision

Reuse the regular mixture transport face fluxes and interpolate them from
Cartesian face centers to EB face centroids. Form the conservative divergence
with open-face areas and each cell's fluid volume. Add no diffusive cut-face
source, defining an adiabatic slip and species-impermeable embedded wall.

Before divergence, compare outgoing species flux with the fluid-volume species
inventory and apply the minimum adjacent-cell positivity factor to the entire
coupled face flux. Advance each transport Euler stage through the existing
StateRedist/EOS transaction, compose two stages as SSPRK2, and insert symmetric
transport half-steps around EB hydrodynamics and inside the chemistry half
steps. Select the public timestep from both hyperbolic and transport limits.

## Consequences

The single-level EB application now shares transport coefficients and flux
physics with the regular solver while retaining EB conservation, small-cell
stability, covered-cell immutability, and transactional rollback. The whole
face scaling is deliberately conservative and keeps species enthalpy and
viscous work coupled when positivity limiting activates.

This is not a no-slip, isothermal, catalytic, or prescribed-flux embedded-wall
model. Coarse/fine diffusive flux matching and reflux, dynamic EB AMR transport,
MPI ownership, and parallel transport remain separate work.
