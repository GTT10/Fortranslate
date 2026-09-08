# Decision 0199: Compose regular-grid molecular transport transactionally in 3D

## Context

Milestone `0.205.0` qualified cell-local chemistry around periodic
three-dimensional general-EOS hydro, but deliberately omitted diffusive
physics. The established 2D transport module already owns the mixture
coefficient evaluation, mixture-averaged species flux, barodiffusion,
correction velocity, and species-enthalpy coupling. Reimplementing those
models would create an avoidable coefficient-parity risk; directly treating
3D as independent planes would omit z faces and the off-diagonal viscous
stress terms.

## Decision

Reuse the established face coefficient and species-flux helpers while adding
a direct periodic x/y/z finite-volume transport operator. Evaluate the full
3-by-3 velocity gradient at each face, apply the Newtonian stress with Stokes'
hypothesis, Fourier conduction, and all-species enthalpy flux, and store one
upper face per cell and direction. Limit outgoing species mass over all six
faces before applying the conservative divergence.

Advance transport with a transactional SSPRK2 candidate and select its stable
step from
`transport_cfl / (Dmax * (1/dx^2 + 1/dy^2 + 1/dz^2))`. Compose the complete
field update as `R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)` inside one outer
transaction. Publish state, recovered temperature, and the minimum species
limiter only after every stage succeeds.

Keep the boundary periodic, serial, single-level, and uniform-grid in this
increment. Qualify x/y/z reduction against the established 2D SSPRK2
operator for elementary and full-H2O2 tables, a genuinely 3D diagonal viscous
wave, independent thermal and species smoothing, and a public full-H2O2
transport/control pair.

## Consequences

The regular 3D path now contains the selected molecular-transport subset and
the same full-physics split ordering used by the qualified lower-dimensional
paths. It conserves periodic Euler and individual species integrals when
chemistry is disabled and rolls back after a valid reaction or transport
prefix if hydro rejects the candidate.

This is not full PeleC transport parity. Soret and Dufour effects, bulk
viscosity beyond Stokes' hypothesis, multicomponent diffusion, nonperiodic
walls, high-order 3D hydro, AMR/reflux, restart, MPI decomposition, EB, LES,
particles/spray, accelerators, and external physical validation remain open.
