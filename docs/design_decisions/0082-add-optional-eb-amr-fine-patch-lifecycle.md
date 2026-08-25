# Decision 0082: add an optional EB AMR fine-patch lifecycle

## Context

Solution-driven EB AMR could move and resize its one fine rectangle, but empty
temperature tags always retained that patch. The application therefore kept
fine storage, subcycling, and fine output even after the resolved feature
disappeared. It also had no qualified transition back from a root-only state.

## Decision

Add separate controls for initial regrid evaluation and removal of an untagged
fine patch. Keep removal disabled by default, and require dynamic regridding
when it is enabled.

Treat removal as a topology transaction. Average the complete fine patch into
a candidate root with reactive EOS recovery, then commit the root and release
fine state, temperature, geometry, and patch metadata. Dispatch subsequent CFL
selection, hydrodynamic advance, and conserved diagnostics through the existing
single-level EB path.

If temperature tags return, build the planned fine geometry and initialize it
by piecewise-constant prolongation from the synchronized root. Publish the new
fine hierarchy only after geometry validation and active-cell EOS recovery.
Write fine output only while the patch is active.

## Consequences

The serial application can now transition conservatively between a root-only
mesh and one strictly internal ratio-aligned fine rectangle. Removing refinement
also removes its storage and subcycling cost, while the default retains the
previous behavior.

The lifecycle still permits at most one fine rectangle. Multiple patches,
deeper levels, EB AMR chemistry and molecular transport, checkpoint/restart,
and MPI ownership remain outside this decision.
