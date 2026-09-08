# Decision 0195: Establish the three-dimensional directional core

## Context

The qualified flow applications stop at two dimensions even though the base
and reacting state layouts already retain three momentum components. Building
a 3D driver before checking the z-coordinate mapping would mix state-layout,
Riemann, mesh, update, and boundary errors in one milestone.

## Decision

Introduce a uniform 3D cell-center API. Treat the existing x-normal Riemann
implementations as canonical and add reversible z-to-x state rotations plus
x-to-z flux rotations. Apply that contract to ideal-gas, passive-multispecies,
and general-EOS reacting states. Require analytic equal-state flux placement,
species closure, and rotation round trips in focused tests.

## Consequences

All currently supported inviscid state layouts now have tested x/y/z normal
flux semantics without duplicating a Riemann solver. A future 3D update can
assemble coordinate fluxes on this foundation. This decision does not qualify
a 3D timestep, boundary condition, transport operator, AMR hierarchy, restart,
MPI decomposition, EB geometry, or application.
