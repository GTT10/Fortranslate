# Decision 0171: Localize multilevel EB conservation closure

## Context

The qualified EB AMR paths measured a conservative residual after reflux and
average-down, then spread it over every active unrefined parent cell. This was
transactional and conservative but allowed a cut or regular coarse/fine
interface to perturb cells arbitrarily far from that interface.

## Decision

Build one replicated recipient mask from direct child rectangles and parent EB
geometry. Mark the active, unrefined parent cells in clipped three-by-three
neighborhoods of coarse cells immediately outside each coarse/fine side. Union
sibling supports, normalize by their total fluid volume, and apply the existing
density, total-energy, and species correction only on that mask. Physical sides
without an exterior parent cell add no support.

Use the same helper in fixed-depth, multipatch, arbitrary-depth, serial, and
sparse-MPI closures. Preserve EOS recovery, composite validation, collective
acceptance, and rollback.

## Consequences

Residual closure is now confined to the coarse/fine neighborhood and no longer
changes unrelated parent cells. Regular and EB-cut child interfaces share the
same deterministic support rule, which is required by the simplified
arbitrary-depth flux-register path. The implementation does not claim exact
AMReX `MLStateRedistribute` per-neighborhood movement bookkeeping.
