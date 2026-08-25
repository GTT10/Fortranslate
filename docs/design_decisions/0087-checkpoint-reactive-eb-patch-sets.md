# Decision 0087: checkpoint reactive EB patch sets separately

## Context

The public reactive EB AMR application could own multiple separated children,
but its established checkpoint schema encoded only one optional fine patch.
Extending that record in place would either break existing readers or require
ambiguous interpretation of its fine-active flag. Restart also has to restore
the complete ordered set atomically because publishing only some children
would invalidate composite ownership and regrid cadence.

## Decision

Keep the single-patch magic and schema unchanged. Add a distinct versioned
formatted patch-set checkpoint whose header declares the child count. Store the
same species, mesh, EB, chemistry, hydro, redistribution, time, and regrid
compatibility information, plus the collection gap setting. Store the complete
root followed by every deterministic child's actual coarse bounds, dimensions,
conserved state, and temperature.

On read, rebuild the root and every child geometry from the compatible input
and stored bounds. Treat conserved state as authoritative and recover active
temperatures through the EOS. Hold all arrays, geometry, metadata, time, and
counters in private candidates until child separation, dimensions, finite
state, EOS recovery, and the terminal marker all validate.

Connect scheduled writes to the public patch-set time loop after the accepted
physics interval and any due periodic regrid commit. Preserve stop-after-write
and final-checkpoint behavior, and resume from stored step and regrid counters.

## Consequences

Serial reacting EB runs can stop and resume with zero, one, or multiple active
children while preserving topology, state, deterministic order, and cadence.
Malformed or incompatible files expose no partial hierarchy. Existing
single-patch checkpoints remain readable with their original schema.

The two formats are intentionally distinct. Cross-format conversion,
distributed EB checkpoint I/O, deeper levels, and EB AMR molecular transport
remain outside this decision.
