# Decision 0198: Compose cell-local chemistry transactionally in 3D

## Context

Milestone `0.204.0` qualified the full reactive state and NASA7 thermodynamic
closure through periodic three-dimensional Euler evolution, but deliberately
advanced no reaction source. The lower-dimensional code already has qualified
elementary explicit and full-H2O2 implicit constant-volume cell reactors and a
transactional chemistry operator. Reimplementing those kinetics inside a 3D
loop would create an unnecessary second numerical path.

## Decision

Advance a private three-dimensional candidate one z plane at a time through
the established transactional 2D cell-local chemistry kernel. Commit the
caller state and temperature only if every plane succeeds. Compose that
operator as `R(dt/2)-H(dt)-R(dt/2)` inside a second whole-step transaction so a
late chemistry or hydro rejection cannot expose a partially reacted field.

Keep chemistry optional and retain the chemistry-disabled entropy-wave path.
Add `uniform_reactor` for exact independent-cell reduction and a periodic
constant-pressure Gaussian `reactive_hotspot` for coupled application testing.
Conserve the first five Euler integrals across a reacting run and verify
species conversion through H/O/N elemental totals rather than incorrectly
requiring individual species integrals to remain fixed.

## Consequences

The serial uniform-grid 3D path now supports reacting elementary and full-H2O2
states with the same cell integrators as the qualified 1D/2D implementations.
The public hotspot is a short-time coupling and conservation regression, not
an ignition-delay or external-validation result. Molecular transport,
arbitrary mechanism ingestion, production stiff-integration parity,
high-order 3D reconstruction, AMR, restart, MPI, EB, LES, particles, and
whole-application PeleC parity remain open.
