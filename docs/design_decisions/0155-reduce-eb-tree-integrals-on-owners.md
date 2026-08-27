# Decision 0155: Reduce EB tree integrals on owners

## Context

Recursive hydrodynamic conservation closure needs a composite conserved
integral before and after every refined subtree advance. The sparse numerical
tree cannot use the serial operation directly because nonowners intentionally
hold no state or temperature fields.

## Decision

Keep topology traversal replicated and deterministic, but integrate a node's
unrefined cells only on its owner. Build the refined mask from direct child
rectangles before visiting descendants. Sum one rank-local conserved vector
over the communicator and separately sum the number of contributing owner
nodes.

Make the full-tree operation a root-subtree wrapper. Require communicator
consensus for the selected level, patch, and vector extent, and publish neutral
outputs on any invalid or rank-dependent selector. Commit optional local-node
accounting only with a finite, nonempty global result.

## Consequences

Whole-tree and arbitrary-subtree conservation measurements no longer
materialize numerical fields. The reduction preserves the serial composite
coverage rule and gives owner-local recursive hydro a collective conservation
boundary. Floating-point summation may differ from serial traversal order, so
parity is qualified within roundoff rather than bitwise equality.
