# Decision 0172: Add embedded-wall diffusive fluxes

## Context

The EB molecular-transport operator interpolated Cartesian diffusive fluxes to
open-face centroids but assigned zero diffusive flux to the cut face. Every
embedded wall was therefore adiabatic, free slip, and species impermeable even
when the existing domain-wall model supported temperature and velocity data.
PeleC similarly reserves nonzero embedded diffusive flux for isothermal heat
transfer and no-slip momentum transfer.

## Decision

Store one embedded-wall record in the existing reactive 2D boundary set. Keep
the default stationary, adiabatic, slip, and impermeable. In each cut cell,
measure the physical fluid-centroid to boundary-centroid distance along the
solid-to-fluid normal. Evaluate mixture viscosity and conductivity at the
recovered cell state, apply a first-order normal Fourier gradient for an
isothermal wall, and apply the normal-only Newtonian stress for a no-slip wall.
Use wall velocity for viscous energy work.

Return a conservative wall-normal flux with exact zero mass and species
components. Multiply it by embedded-wall length and divide by cut-cell fluid
volume in the existing transport right-hand side. Reuse the unchanged
StateRedist, AMR, sparse-MPI, EOS, consensus, and rollback paths.

## Consequences

Every library path that already receives a validated boundary set can use the
same thermal/no-slip wall flux without a new AMR or MPI interface. MPI control
matching now includes the embedded record. The model is first order and does
not reproduce PeleC's quadratic EB boundary-gradient stencil. Catalytic or
prescribed species transfer, namelist exposure, and checkpoint fingerprints
for nondefault embedded walls remain later work.

