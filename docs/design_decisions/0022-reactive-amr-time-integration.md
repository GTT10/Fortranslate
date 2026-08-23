# Decision 0022: first runnable reactive AMR time integration

## Status

Accepted for PeleF `0.27.0`.

## Context

The hierarchy, transfers, flux register, tagging, and regrid planner were
individually verified but did not yet advance a fluid calculation. A useful AMR
milestone must execute a reacting case and apply the synchronization primitives
inside its time loop.

## Decision

Introduce a dedicated two-level reactive solution and executable. Advance PCM
general-EOS Godunov hydro once on the coarse level and refinement-ratio times on
the fine patch. Interpolate adjacent coarse conserved states in time for fine
ghost cells. Accumulate both interface fluxes over the coarse and fine time
steps, reflux, and average down.

Apply cell-local chemistry for half an interval on both levels before hydro and
for half an interval after synchronization. Average down again after the final
chemistry operation. Make the full hierarchy step transactional and evaluate
the existing solution-driven regrid planner at a configured coarse-step
interval.

## Consequences

- AMR now advances an actual reacting flow rather than isolated transfer data;
- composite mass, momentum, and energy have executable conservation gates;
- output contains a non-overlapping ordered composite mesh;
- PCM limits the first integration slice to a robust flux-consistency baseline;
- AMR molecular transport, high-order coarse/fine reconstruction, multiple
  patches, multilevel recursion, and MPI ownership remain follow-on work.
