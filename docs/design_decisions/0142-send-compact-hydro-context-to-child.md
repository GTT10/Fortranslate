# Decision 0142: Send compact hydro context to each child

## Context

Sparse hydro advances root bands on tile owners and assembles their result on
the root physics owner. That owner previously sent complete root start/end
state, temperature, and x/y flux arrays to every distinct remote child owner.
Each child needed only its exterior interpolation edges, nearby reflux state,
and coarse/fine interface fluxes.

## Decision

For each child, extract the established four-edge stage-start/stage-end
context, current patch-plus-two corrected state and temperature, and the
patch-boundary x/y flux rectangle on the root physics owner. Pack them into one
message to a remote child. Build fine exterior state from the compact context,
accumulate the coarse register through the globally indexed flux-support API,
and reflux through the globally indexed state-support API on the child owner.
Return only corrected patch-plus-two state and temperature to the root owner.

Preserve deterministic child order and merge each corrected support before the
next child. Keep root-tile result assembly and final corrected-row scatter in
this milestone. Hydro and transport retain separate message tags.

## Consequences

Remote hydro child owners allocate no complete root state, temperature, or
flux field. A remote child transaction contains two messages regardless of the
complete root size. Qualification must require a strictly smaller payload,
exact traffic, complete/support numerical parity, owner accounting,
conservation, and rollback in Debug and Release at one, two, four, and eight
ranks. This decision does not claim a measured speedup.
