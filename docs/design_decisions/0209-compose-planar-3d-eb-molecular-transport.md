# Decision 0209: Compose molecular transport on the planar 3D EB state

## Context

The regular three-dimensional solver already evaluates Newtonian viscosity,
Fourier conduction, mixture-averaged species diffusion, barodiffusion,
correction velocity, and species-enthalpy flux. The two-dimensional EB path
also establishes the required aperture weighting, species-outflow limiting,
wall transport, and StateRedist composition. The `0.215.0` planar 3D EB
application has neither a transport timestep nor a transport operator, so its
`R-H-R` split cannot represent the selected reacting-flow physics set.

The first 3D EB geometry is deliberately narrower than the 2D geometry. It is
one non-face-aligned x-, y-, or z-normal plane, with a single regular receiver
for every small cut cell. Extending transport at this boundary is useful only
if it preserves that restriction instead of implying general cut-face
gradients or wall models.

## Decision

Add a dedicated planar 3D EB transport module. Recover all active-cell
primitives transactionally, set every physical-domain outer transport flux to
zero, evaluate the existing multicomponent constitutive laws only on open
interior Cartesian faces, and weight their divergence by stored apertures and
physical fluid volume. The embedded wall is fixed to adiabatic slip and
impermeable, so its viscous, conductive, and species transport flux is also
exactly zero.

Limit outgoing species flux with each cell's physical fluid inventory before
forming the divergence. Advance each transport interval with transactional
SSPRK2, applying the qualified zeroth-order StateRedist update at both Euler
stages. Covered state and temperature storage must remain bitwise identical.
Use an active-state full-domain diffusive bound, without a separate `kappa`
penalty because both Euler stages use StateRedist, and select every public step
as the minimum of hydro and transport limits.

Extend the EB split to `R-T-H-T-R`: reaction and transport receive half of the
accepted interval on either side of one stabilized hydro step. All stages
remain private until the complete candidate passes. With chemistry or
transport disabled, skip the corresponding stages without perturbing the
existing arithmetic. Public transport requires StateRedist; FluxRedist
transport remains unqualified.

## Consequences

This increment supports the existing elementary mixture on one serial,
single-level, exact-axis-plane mesh with PCM hydro, zero-gradient outer faces,
and an adiabatic slip/impermeable embedded wall. It can test viscous,
conductive, and mixture-averaged species transport without introducing a new
wall heat or mass source.

The deterministic public control, transport-only, and coupled cases stop after
one clipped `2.0e-9 s` step, before the zero-gradient hydro boundary can become
part of the conservation comparison. They qualify the complete operator
transaction and a resolved local transport/chemistry response; they are not a
long-time boundary-flow validation.

Public invariant diagnostics use the maximum initial/final L1 inventory for
each component. This keeps the reported relative error meaningful when the
same problem is scaled to a small physical volume; it does not excuse genuine
exchange once a disturbance reaches the zero-gradient hydro boundary.

It does not qualify isothermal or no-slip embedded walls, prescribed species
flux, catalytic chemistry, general geometry, higher-order face-centroid
interpolation, FluxRedist transport, AMR/reflux/regridding, restart, MPI
ownership, configurable mechanisms, or external physical validation. The
step remains selected before the first reaction half step, so the coupled
public case is a frozen low-feedback regression rather than a general
reaction-aware stability result.
