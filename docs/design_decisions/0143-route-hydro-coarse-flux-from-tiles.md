# Decision 0143: Route hydro coarse flux from tile owners

## Context

Sparse hydro advances finite root bands on tile owners. The root physics owner
still needs complete stage-start and stage-end state and temperature for child
exterior extraction, ordered reflux-support merge, and final row scatter. It
does not need a complete x/y flux field once child registers can consume
globally indexed interface rectangles.

The compact hydro child transaction introduced in `0.150.0` still moved every
child's coarse interface flux through the root physics owner. Remote tile
results also carried owned flux rows solely to build that temporary root flux
bundle.

## Decision

Retain each hydro tile's x-flux rows and uniquely owned y-faces on the tile
owner. For each child, route one packed x/y fragment from every intersecting
remote tile owner directly to the child owner. The receiver assembles the
globally indexed interface rectangles and requires complete finite coverage
before accumulating the coarse register.

Remove flux from the tile-to-root hydro result and from the root-to-child state
context. Keep the existing state-only root assembly, child-local subcycling and
reflux, ordered corrected-support return, and final corrected-row scatter.
Preserve the unique y-face rule so adjacent tiles neither omit nor duplicate an
interface face.

## Consequences

The root physics owner allocates no complete hydro x/y flux array. Hydro flux
traffic scales with tile/child boundary intersections rather than the complete
root face count, and a child-local tile needs no message. Qualification derives
the exact message count from those intersections plus existing halo, state
context, corrected-support, tile-result, and final-scatter routes. Serial field
parity, conservation, rollback, and owner accounting remain required in Debug
and Release at one, two, four, and eight ranks. This decision does not claim a
measured speedup.
