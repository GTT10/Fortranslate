# Decision 0132: advance sparse transport root tiles with finite halos

## Context

After the tile-local final SSPRK2 blend, both sparse root transport Euler stages
still gather every root tile to one physics owner before any flux or
redistribution work begins. Diffusive face construction, flux limiting, EB face-
centroid interpolation, and second-order StateRedist need neighboring values,
so the root update cannot be split into independent zero-halo tiles.

The established owner-tiled hydro path demonstrates the required ownership
boundary: assemble a bounded geometry/state band on each target owner, execute
the serial kernel unchanged on that band, and assign deterministic ownership to
the result and shared flux rows. Transport has a composed dependency from
face fluxes through regular and EB limiting, the conservative RHS, overlapping
redistribution neighborhoods, and second-order neighborhood-state slopes. The
established hydro decomposition already qualifies a conservative six-row guard
for the face-centroid and second-order StateRedist composition.

## Decision

For each sparse SSPRK2 Euler stage, give every root tile a six-row guard.
Copy same-owner fragments directly and send only intersecting remote row
fragments point to point. On the target owner, execute the qualified EB
transport-flux RHS and StateRedist kernels on the complete band. Route the
tile's owned input, updated state, recovered temperature, x-flux rows, and
uniquely owned y-faces to the root physics owner.

If the target touches a periodic y boundary, assemble the complete root band.
The serial transport kernel wraps the ends of its supplied geometry, so a
clipped edge band would connect the physical boundary to the wrong interior
row. This compatibility case remains explicit until cyclic band geometry is
available.

Retain the complete temporary root bundle after tile computation for the
existing fine-child exterior, diffusive flux register, and cumulative reflux
boundary. Scatter only the corrected final row bands. Count halo, result,
child-bundle, correction, and scatter messages exactly, and separately count
every band cell passed through the transport kernel. Publish neither count
until the complete SSPRK2 transaction succeeds.

## Consequences

No unconditional gather to one selected owner occurs before either Euler-stage
transport advance. Persistent fields remain globally single-copy, and interior
targets use work determined by their tile plus a fixed guard. Small roots and
periodic boundary targets can use a complete temporary band; cyclic finite-halo
geometry is needed before claiming a fixed-width boundary for those targets.

The root physics owner still receives one post-compute level-wide bundle per
stage. Removing that allocation and the child correction round trips requires
distributed coarse/fine exterior and reflux data, which remains separate work.
