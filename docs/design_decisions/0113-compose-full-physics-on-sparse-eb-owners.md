# Decision 0113: compose full physics on sparse EB owners

## Context

The sparse EB full-physics entrypoint retains owner-only input and output but
still materializes a complete hierarchy for its central
transport-hydro-transport window. Direct sparse chemistry, hydro, and SSPRK2
transport are now qualified independently, so that compatibility window no
longer protects an unsupported component.

## Decision

Compose one private sparse candidate as direct chemistry, direct SSPRK2
transport, direct hydro, direct SSPRK2 transport, and direct chemistry. Pass
the same immutable geometry and patch-set template to each operator, and retain
the caller's sparse object until all five component transactions succeed.

Publish chemistry, hydro, and transport counts and the minimum transport
limiter only with the final sparse commit. If any later component rejects,
discard the candidate and retain the caller's fields and fallback bookkeeping
exactly.

## Consequences

The complete `R-T-H-T-R` split no longer allocates or synchronizes a complete
fine-child hierarchy. Fine numerical payloads stay globally single-copy across
operator boundaries, while each component preserves its independently
qualified root-level temporary and collective acceptance rules.

Root physics decomposition, targeted root and coarse/fine communication,
public time-loop control, dynamic topology, checkpointing, and output remain
future work.
