# Decision 0207: expose a stabilized planar 3D EB application

## Context

The 3D EB work through `0.213.0` was callable only from library regressions.
Its stabilized algorithms intentionally accept a full Cartesian CFL step, but
there was no public timestep selector, input contract, repeated-step driver,
or field output carrying the EB metrics needed for independent audit.

A public raw cut-cell route would make the default application unstable by
construction. Adding chemistry, transport, AMR, or general geometry at the
same boundary would also obscure whether failures came from stabilization or
from newly coupled physics.

## Decision

Add `compute_reactive_eb_cfl_timestep_3d`. It skips covered cells, recovers
every active NASA7 primitive state, and selects

```text
dt = CFL / max[(|u|+c)/dx + (|v|+c)/dy + (|w|+c)/dz].
```

Do not divide this full-grid rate by cut-cell volume fraction; the public
driver requires either `state_redist` or `flux_redist`. Reject an all-covered
domain, invalid active state, bad extents, or CFL outside `(0,1]`.

Add a dedicated `&reactive_eb_3d` input record and the installed
`pelef_reactive_eb_3d` executable. The configuration admits only one exact
axis-normal plane with at least one regular receiver, PCM Rusanov/HLLC/PeleC
fluxes, and the two stabilized routes. Public StateRedist additionally
requires its target to exceed the constructed cut fraction, so the named
merge cannot silently become an unstabilized update. The raw library routine
remains a test control and cannot be selected publicly.

Initialize a frozen elementary H2/O2/N2 mixture at constant pressure with one
density on regular/covered cells and another on the cut sheet. Advance
transactional forward-Euler candidates until `final_time`, clipping the final
step. Publish x-fastest CSV containing Cartesian indices and centers, cell
type, volume fraction, physical fluid centroid, primitive state, and all
species conserved fields.

Treat wall-normal momentum as a physical exchange with the embedded wall, not
as a conserved fluid invariant. Gate mass, total energy, all species, and the
two tangential momenta. Report the wall-normal momentum change separately.

## Consequences

The default `10 x 8 x 6`, `kappa=0.05`, CFL-0.5 case takes two accepted steps
to `5.0e-5`. StateRedist and FluxRedist conserve their fluid invariants to
roundoff but produce intentionally distinct density and wall-impulse
signatures. An independent Python checker reconstructs all 480 geometry rows,
EOS/composition relations, fluid-volume integrals, resolved response, and
cross-method invariant agreement.

This is the first public 3D EB hydro application, not a production PeleC
replacement. It is serial, single-level, inviscid, first-order in time and
space, zero-gradient at the Cartesian boundary, elementary-mixture only, and
exact-axis-plane only. Chemistry, transport, characteristic reconstruction,
general geometry, AMR/reflux/regridding, restart, MPI ownership, scalable I/O,
and external physical validation remain separate work.
