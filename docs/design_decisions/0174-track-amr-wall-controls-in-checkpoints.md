# Decision 0174: Track AMR wall controls in every restart contract

## Context

The `0.181.0` public namelist exposed embedded-wall controls only to the
checkpoint-free single-level application. Public AMR applications rejected an
active isothermal or no-slip selection because their fixed-depth checkpoint
schemas and arbitrary-depth fingerprint did not identify the wall or the
molecular-transport controls that make its flux meaningful.

## Decision

Build the boundary set through one configured helper in every public EB driver.
Advance the fixed-depth two-level, sibling-patch, static three-level, and dynamic
three-level checkpoint schemas to version 2. Store and compare wall kind,
thermal mode, temperature, velocity, transport enable, all transport process
switches, and transport CFL.

Advance the arbitrary-depth checkpoint fingerprint to schema 3 and add the same
wall and transport identity. Share that fingerprint between serial and
sparse-MPI checkpoint paths. Reject older schemas or any changed compatibility
value transactionally instead of substituting defaults.

## Consequences

Checkpoint-free and restartable AMR applications can now activate the same
embedded heat and momentum transfer as the single-level application. A restart
cannot silently change wall physics or transport selection. The schema changes
are intentionally incompatible with older formatted checkpoints; state-file
migration is not part of this milestone.
