# Decision 0068: use a pressure-only reactive EB slip-wall flux

## Context

The embedded-boundary geometry stores a unit normal from solid toward fluid,
while a finite-volume surface flux uses the outward normal of the fluid control
volume. A stationary inviscid impermeable wall transfers pressure force but no
mass, energy, or chemical species through the wall.

## Decision

Recover pressure from the full reactive conserved state with the qualified
general-EOS conversion. For stored fluid normal `n`, return an outward wall
flux whose only nonzero entries are the in-plane momenta `-p*n`. Integrate that
flux over the physical cut-interface length and subtract it from the control
volume balance, giving a source `p*n*A/(kappa*dx*dy)` in each cut cell.

Require a finite unit normal, positive temperature guess, valid reactive state,
positive cut volume, and positive interface length. Build the complete source
in temporary storage so any rejected cut cell leaves the public output zero.

## Consequences

The kernel is rotationally covariant and conserves the impermeability contract
for arbitrary planar wall orientation. It does not yet combine Cartesian
open-face fluxes, advance a cut-cell state, stabilize small cells, or model
viscous, thermal, moving, or catalytic wall effects.
