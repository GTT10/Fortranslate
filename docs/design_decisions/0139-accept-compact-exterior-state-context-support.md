# Decision 0139: accept compact exterior state-context support

## Context

Direct coarse-flux routing removed flux-register construction from the root
physics owner, but child exterior-context extraction still required complete
root start/end state and temperature arrays. The extractor reads only coarse
cells immediately outside the four patch sides.

## Decision

Add an exterior-context extraction entrypoint whose coarse arrays carry
explicit global lower bounds. Require a rectangle containing the fine patch
expanded by one coarse cell, clipped to the root domain. Validate bounds,
shapes, and finite values before publishing the context.

Keep the complete-root extraction API as a wrapper over the support kernel.
Require bitwise-identical reconstructed edge state and temperature from a
strictly smaller support in the unit gate.

## Consequences

Context extraction no longer requires complete-root array shapes. Sparse MPI
still uses the root wrapper until start/end state fragments and cumulative
corrected support move between root tile and child owners. Interpolation,
physical-side placeholders, EOS recovery, and child subcycling are unchanged.
