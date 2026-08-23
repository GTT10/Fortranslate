# Decision 0020: static conservative AMR 1D foundation

## Status

Accepted for PeleF `0.25.0`.

## Context

PeleF needs AMR transfer and synchronization rules before reactive operators or
MPI ownership can be attached to refined patches. Reimplementing the complete
AMReX hierarchy would obscure the smaller subset required by the porting plan.

## Decision

Begin with one coarse level and one strictly interior refined patch. Store level
and index-box metadata independently of the state width. Require an integer
refinement ratio of at least two.

Use MC-limited piecewise-linear conservative prolongation, volume-average
restriction, refinement-ratio subcycling, and a two-interface flux register.
Accumulate time-integrated fine-minus-coarse interface fluxes and reflux only
the adjacent uncovered coarse cells. Count uncovered coarse volumes and refined
volumes exactly once in composite diagnostics.

## Consequences

- the transfer routines work for Euler, multispecies, or reactive state widths;
- prolongation followed by restriction preserves coarse averages;
- reflux has a direct composite-conservation proof and regression;
- fine patches touching physical/periodic domain boundaries are deferred;
- tagging, regridding, more than two levels, solver coupling, and MPI patch
  ownership remain explicit follow-on decisions.
