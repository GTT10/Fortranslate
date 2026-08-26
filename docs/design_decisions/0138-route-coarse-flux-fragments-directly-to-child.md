# Decision 0138: route coarse flux fragments directly to child

## Context

Compact coarse flux-register accumulation removed the numerical requirement
for complete root flux arrays, but sparse transport still initialized every
child register on the root physics owner and sent it with the state context.
Root tile owners already held all required interface values immediately after
their bounded transport computation.

## Decision

Retain each root tile's computed x-flux rows and uniquely owned y-faces with
global y bounds. For each child, send one packed fragment from every remote
intersecting tile owner directly to the child owner. Assemble the complete
compact x/y face rectangles there and reject collectively unless every face is
covered by the deterministic ownership partition.

Initialize and accumulate the coarse register on the child owner. Remove the
register from the root-to-child exterior/state-support context. Preserve
child-local fine accumulation, reflux, corrected-support return, cumulative
root merge order, and transactional publication.

## Consequences

The root physics owner no longer creates or communicates coarse child
registers. Message count now depends on root-tile/child intersections, while
each payload is restricted to actual interface fragments. The root still
assembles complete temporary result and flux bundles for state context,
boundary closure, support merge, and final scatter; those remain later
distribution boundaries.
