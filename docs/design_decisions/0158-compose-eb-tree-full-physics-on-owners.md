# Decision 0158: Compose EB tree full physics on owners

## Context

Sparse chemistry, recursive hydro, and recursive SSPRK2 transport now each
commit atomically on numerical-node owners. Calling them directly on the
accepted tree would still expose a valid reaction or transport prefix if a
later split stage failed.

## Decision

Own one outer sparse candidate and apply optional chemistry, transport, hydro,
transport, and optional chemistry in `R-T-H-T-R` order. Make timestep,
tolerances, redistribution controls, physics flags, and data extents
communicator-consistent before optional branches. Continue using each inner
operator's boundary and scheme consensus.

Keep chemistry, transport, and hydro advances and transfer counts in separate
private categories. Combine the two transport limiter minima. Assign the
candidate and publish all diagnostics only after every stage and final sparse
validation succeed.

## Consequences

The normal arbitrary-depth full-physics path never materializes a complete
tree and cannot publish a partial split prefix. Its traffic is exactly two
chemistry restriction schedules, two SSPRK2 transport schedules, and one hydro
schedule when every operator is enabled. A target-time loop can now use this
operation as its per-step transaction without weakening rollback.
