# Decision 0090: synchronize three-level EB hierarchies deepest first

## Context

The qualified EB AMR transfer operated on one root and one fine rectangle.
Advancing toward deeper EB AMR first requires an ownership rule for composite
diagnostics and a synchronization order that cannot overwrite finest-level
information before it reaches the root.

## Decision

Represent the initial deeper-level foundation as two strictly nested aligned
patch descriptors: root to middle and middle to finest. Count root cells only
outside the middle patch, middle cells only outside the finest patch, and all
finest cells in composite integrals.

Synchronize deepest first. Average the finest state into a private middle
candidate, then average that complete candidate into a private root candidate.
For reactive state, recover an EOS-consistent temperature after each parent
restriction. Validate all arrays and geometry up front and publish both parent
levels only when both stages succeed.

## Consequences

Finest data reaches the root without double counting or order-dependent loss,
and a failure leaves both parent levels unchanged. The implementation reuses
the established cut-cell volume weighting and reactive EOS transaction rather
than introducing a second restriction formula.

This decision does not add hierarchy allocation, prolongation, ghost filling,
recursive hydro subcycling, a second flux register, chemistry composition,
regridding, checkpointing, output, transport, arbitrary depth, or MPI
ownership. Those remain separate milestones built on this synchronization
contract.
