# Decision 0197: Qualify general-EOS reactive Euler in 3D before chemistry

## Context

Milestone `0.203.0` qualified periodic three-dimensional conservation and
SSPRK2 integration only for the constant-`gamma` Euler state. Extending that
path with reacting thermodynamics, cell-local chemistry, molecular transport,
and AMR in one change would make a failed temperature recovery, species
closure error, source-integration error, or directional divergence error
difficult to distinguish.

## Decision

Add a serial, single-level, uniform-grid, periodic general-EOS Euler path over
the established reactive state layout. Reuse the qualified x/y/z
Rusanov/HLLC/PeleC-style fluxes, store upper-face fluxes directly, and use the
same summed directional CFL and transactional SSPRK2 composition as the
constant-`gamma` path. Recover every candidate cell's NASA7 temperature before
publishing a stage. Reject the complete update without changing either state
or temperature if any face or cell is inadmissible.

Keep chemistry and molecular transport disabled in this increment. Qualify
all Euler and species conserved components, not only the first five. Use an
independent one-dimensional SSPRK2 reference assembled from the public
general-EOS x flux because the older 1D PCM application intentionally uses a
different Forward-Euler time boundary.

Expose `pelef_reactive_3d` with elementary or full-H2O2 thermodynamic tables, a
constant-composition diagonal entropy wave, and deterministic x-fastest CSV
containing both `Y_k` and `rhoY_k`.

## Consequences

The regular-grid 3D path now evolves a composition-dependent NASA7 mixture and
conserves every species density to roundoff. This is a reactive-state
hydrodynamic qualification, not a reacting-flow claim: no reaction source is
advanced. Chemistry splitting is the next isolated increment, followed by
molecular transport. High-order 3D reconstruction, AMR, restart, MPI
decomposition, EB, LES, particles, and production I/O remain open.
