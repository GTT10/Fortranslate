# Decision 0183: Schedule and restart fixed three-level parent regrids

## Context

The atomic parent-patch rebuild exists as a library operation, but the public
fixed three-level application still schedules only its finest patch. Enabling
the parent operation implicitly for every existing dynamic input would change
qualified topologies and make schema-3 dynamic checkpoints ambiguous because
they identify only a fixed configured parent.

## Decision

Add an opt-in `dynamic_parent_regridding` control. It is valid only with
dynamic fixed three-level AMR and with a minimum tagged parent size large
enough to contain an eight-cell middle level plus the finest two-cell margin.
At initialization and each accepted regrid interval, attempt the parent
transaction first. If the parent is unchanged, retain the established
finest-only regrid; if it changes, its transaction already rebuilds the
finest patch and counts as one regrid event.

Advance the dynamic three-level checkpoint to schema 4. Store and match the
new policy flag, accept the actual root-to-middle descriptor when the flag is
enabled, validate the finest descriptor against that stored parent size, and
rebuild both refined geometries privately before publishing restart state.

## Consequences

Existing dynamic three-level inputs keep their fixed parent by default. New
inputs can move or resize both refined levels on the public schedule and
restart that exact topology. A policy mismatch or malformed stored nesting is
rejected transactionally. The fixed-depth path still owns one parent and one
finest rectangle; sibling parent patches require the multipatch or
arbitrary-depth topology.
