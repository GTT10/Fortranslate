# Decision 0016: extend the qualified mixture transport by direction

## Context

PeleF 0.16.0 qualified dilute-gas viscosity, thermal conduction, and
mixture-averaged species diffusion in one dimension.  The reactive two-
dimensional solver already used a full-conserved-state CTU hydro update.

## Decision

Evaluate transport coefficients at cell faces from arithmetic face states,
build Newtonian stress, Fourier, species, and species-enthalpy fluxes in both
coordinate directions, and advance their conservative divergence with SSPRK2.
A single face coefficient limits all species diffusion fluxes when trace-
species positivity would otherwise be lost, preserving zero net diffusive mass
flux and the matching enthalpy transport.  Compose reaction, transport, and
hydro symmetrically.

## Qualification boundary

This is a periodic, uniform-grid, dilute ideal-gas, mixture-averaged subset.
Soret/Dufour effects, bulk viscosity, multicomponent diffusion, physical wall
fluxes, AMR, EB, and parallel execution remain outside the milestone.
