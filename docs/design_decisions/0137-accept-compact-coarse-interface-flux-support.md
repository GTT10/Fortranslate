# Decision 0137: accept compact coarse interface-flux support

## Context

Sparse root transport computes bounded owner-tile bands, but the coarse EB
flux-register API required complete root x/y flux arrays. The root physics
owner therefore remained the only place that could initialize a child's
register, even though accumulation reads only the four coarse/fine interface
face ranges.

## Decision

Add a coarse accumulation entrypoint whose x-face and y-face arrays carry
explicit global lower bounds. Require the supplied rectangles to contain every
active vertical and horizontal coarse/fine interface face, and validate all
bounds and values transactionally. Keep the complete-root API as a wrapper.

Use the compact entrypoint in sparse MPI transport, passing only each child's
interface rectangles from the temporary root flux arrays. Retain the existing
register representation, correction arithmetic, child ordering, and message
route.

## Consequences

The accumulator consumes data proportional to child interface area rather
than root area and produces a bitwise-identical correction. Incomplete or
nonfinite support cannot partially modify a register. The root physics owner
still assembles the complete temporary flux bundle, but the accumulation API
no longer prevents later direct root-tile-to-child interface routing.
