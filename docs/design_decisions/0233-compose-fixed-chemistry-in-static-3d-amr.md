# 0233: Compose fixed chemistry in static 3D AMR

## Status

Accepted for `0.241.0`.

## Context

The `0.240.0` serial and rank-local sparse MPI paths qualified a static,
strictly interior, periodic two-level 3D AMR hierarchy for PCM and
characteristic-PLM hydrodynamics. The regular 3D solver already had
transactional cell-local elementary and full-H2/O2 chemistry, but applying it
independently to two AMR levels would leave covered coarse species and
temperature inconsistent with the fine solution. The sparse MPI path also
needs reaction and integrator metadata consensus before any rank enters a
data-dependent source solve.

## Decision

Add hierarchy-wide fixed-mechanism chemistry and a transactional Strang split
to the existing static 3D AMR boundary.

- A chemistry phase advances private coarse and fine candidates cell by cell,
  averages fine conserved state onto covered coarse cells, and recovers the
  covered coarse temperature before publishing either level.
- A complete reactive step is `R(dt/2)-H(dt)-R(dt/2)` over one private
  hierarchy. Failure of either source phase, hydro, reflux, average-down, or
  temperature recovery leaves caller state and the reported reflux unchanged.
- The sparse MPI source owns only rank-local coarse slabs and parent-aligned
  fine planes. Ranks with no fine cells skip the fine source solve but remain
  in every collective acceptance decision.
- MPI operations compare the complete ordered reaction records, fixed NASA7
  data, patch/ownership, tolerances, integrator policy, solver,
  reconstruction, limiter, timestep, and chemistry-enable flag before
  candidate publication. Rank-dependent valid metadata is rejected.
- Chemistry-disabled calls dispatch through the qualified hydro path and keep
  its exact serial/MPI arithmetic.
- Public diagnostics treat density, three momenta, and total energy as the
  conserved Euler subset. Individual species must change in the reacting
  case; H/O/N composite totals are conserved independently.
- Chemistry-enabled checkpoint/restart is rejected. Existing schema 2 does
  not contain the ordered reaction set, source integrator, or tolerances and
  therefore cannot safely resume this operation.

The supported mechanism boundary is the committed elementary or full-H2/O2
model. Configure-time selected mechanisms, molecular transport, physical
coarse boundaries, dynamic topology, and 3D EB AMR are not enabled by this
decision.

## Consequences

The public full-H2/O2 hotspot changes species and temperature while preserving
Euler and H/O/N composite quantities to roundoff. Serial and sparse MPI
outputs are byte-identical at one, two, four, and eight ranks. Direct tests
also cover invalid reactions and integrators, full-split rollback,
rank-dependent reaction and chemistry-enable policies, collective failure on
a rank with no fine cells, and exact chemistry-disabled hydro dispatch.

The hierarchy remains static and periodic, root gather/scatter remains
compatibility I/O, and no chemical checkpoint format or physical PeleC field
parity is claimed.
