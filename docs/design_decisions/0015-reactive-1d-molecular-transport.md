# Decision 0015: qualify molecular transport in 1D before 2D coupling

## Context

PeleF 0.15.0 had inviscid reactive flow with NASA7 thermodynamics and
characteristic PPM/CTU hydro, but no viscosity, heat conduction, or species
diffusion. Directly adding the full PelePhysics transport stack would combine
coefficient generation, diffusion physics, multidimensional gradients, and CTU
coupling in one untraceable change.

## Decision

Introduce a separately qualified one-dimensional ideal-gas transport slice.
Use pinned Lennard--Jones records and standard dilute-gas formulas for pure and
mixture coefficients. Assemble conservative face fluxes in the same sign and
state layout as PeleC `Diffterm.H`:

- Newtonian shear stress and viscous work;
- Fourier conduction;
- mixture-averaged mole-fraction diffusion;
- optional ideal-gas barodiffusion;
- correction velocity enforcing zero net diffusive mass flux;
- species-enthalpy diffusion in the total-energy flux.

Advance this operator with SSPRK2 under an explicit parabolic timestep limit.
Compose chemistry, transport, and hydro symmetrically. Keep transport disabled
by default so every earlier inviscid regression is unchanged.

## Qualification boundary

The coefficient model is not full PelePhysics `SimpleTransport` parity. It
omits generated polynomial fits, polar corrections, bulk viscosity, Soret,
Dufour, and multicomponent Stefan--Maxwell diffusion. It is currently coupled
only to the periodic/outflow one-dimensional reactive path.

## Consequences

Coefficient algebra, conservative diffusion fluxes, analytical convergence,
and reaction-flow coupling can now be tested independently before any 2D
gradient or corner-coupling work. The next transport milestone can reuse these
qualified normal fluxes in the reactive 2D solver.
