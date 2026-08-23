# Decision 0024: conservative molecular transport on reactive AMR

## Status

Accepted for PeleF `0.29.0`.

## Context

The dynamic two-level reactive hierarchy already synchronized advective fluxes
and cell-local chemistry, but rejecting molecular transport prevented viscous,
conductive, and diffusive reacting calculations. A transport integration must
respect the parabolic fine-grid stability limit and conserve every component
across coarse/fine interfaces.

## Decision

Reuse the qualified one-dimensional diffusive face kernel on both AMR levels.
Advance each transport interval with SSPRK2 and return the arithmetic mean of
its two stage fluxes. Advance the fine level with `r^2` substeps, interpolate
coarse ghost states to each substep midpoint, and evaluate interface gradients
over the actual coarse-center to fine-center distance.

Accumulate the stage-averaged coarse and fine diffusive fluxes in a flux
register, then reflux and average down. Compose transport symmetrically with
reaction and hydro as `R(dt/2)-T(dt/2)-A(dt)-T(dt/2)-R(dt/2)`. Keep the complete
hierarchy update transactional.

## Consequences

- viscosity, conduction, and mixture-averaged species diffusion now operate on
  the dynamic two-level reactive hierarchy;
- parabolic subcycling is explicit and scales quadratically with refinement;
- diffusive synchronization conserves every periodic state component and each
  species mass;
- the existing dilute-gas transport qualification boundary is unchanged;
- Soret, Dufour, multicomponent Stefan--Maxwell diffusion, multiple patches,
  arbitrary multilevel recursion, and MPI patch distribution remain follow-on
  work.
