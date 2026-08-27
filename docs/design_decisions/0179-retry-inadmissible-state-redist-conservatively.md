# Decision 0179: Retry inadmissible StateRedist conservatively

## Context

Cut-parent limited-linear prolongation preserves physically admissible child
states and the parent average, but it also retains fine-scale temperature
gradients that PCM removed. In the dynamic three-level reactive regression,
the second SSPRK2 transport Euler stage exposed a separate robustness gap:
component-limited order-2 StateRedist produced a conserved-variable
combination that the EOS could not recover.

Componentwise monotonicity does not imply thermodynamic admissibility because
density, momentum, total energy, and species densities are coupled by the EOS.
Committing part of a redistributed result or clipping one component would also
break the existing whole-level transaction or conservation contract.

## Decision

Build and EOS-check the requested StateRedist candidate transactionally. If an
order-2 candidate is inadmissible, discard the complete candidate and
redistribute the unchanged provisional state again with the established
order-0 weighted kernel. EOS-check every active retry cell before committing
either state or temperature.

Both orders use the same overlapping neighborhoods and partition weights, so
each candidate independently preserves every volume-weighted conserved
component. An inadmissible order-0 candidate remains a hard failure and leaves
the caller's output unchanged. Order 0 never retries itself.

## Consequences

Smooth admissible order-2 updates are unchanged. Reactive EB transport and
hydro callers gain a conservative positivity fallback when higher-order
StateRedist combines otherwise bounded components into an inadmissible state.
The fallback is intentionally whole-level and may add diffusion only on a
step that could not otherwise commit. This does not claim a convex
thermodynamic limiter for individual StateRedist neighborhoods.
