# Decision 0084: checkpoint the reactive EB AMR lifecycle

## Context

The serial reactive EB AMR application could move, remove, and re-create its
fine patch and advance chemistry, but every invocation initialized from the
problem definition. A stopped calculation lost its evolved state, actual patch
bounds, root-only lifecycle state, and accepted-step/regrid cadence.

## Decision

Add a versioned formatted checkpoint with a magic header and terminal marker.
Store species ordering, a strict physics/topology compatibility signature,
actual patch bounds, time and counters, minimum accepted timestep, base density,
and the complete root and optional fine state and temperature payloads.

Write scheduled checkpoints only after the accepted step and any regrid
transaction commit. Permit an input to stop immediately after that successful
write. On restart, rebuild EB geometry from the current input and stored patch
bounds, treat conserved state as authoritative, and recover active temperatures
through the general EOS. Read into private candidates and publish nothing until
the entire file, geometry, EOS state, and end marker are valid.

Keep final time, maximum step count, output paths, and checkpoint controls
outside the compatibility signature so a continuation can extend the run or
change its output schedule. Require mesh, EB geometry, refinement, chemistry,
hydro, redistribution, and dynamic-regrid settings to match.

## Consequences

An active two-level or root-only serial EB AMR run can now stop and resume while
preserving topology and cadence. Malformed or incompatible input cannot expose
a partially restored hierarchy.

The schema is formatted and serial. Distributed EB checkpoint I/O, multiple
fine patches, deeper levels, and EB AMR molecular transport remain outside this
decision.
