# Decision 0080: add a runnable static EB AMR driver

## Context

The two-level reactive EB hydrodynamic API can advance one interval, but callers
still had to build both geometries, construct the patch, initialize the fine
state, select a stable hierarchy timestep, repeat the update, and write output.
That made the AMR path a tested kernel rather than a usable application.

## Decision

Add an `eb_amr` namelist containing a strictly internal coarse rectangle, an
integer refinement ratio, and a distinct fine output path. Reuse the existing
reactive EB configuration for the root grid, geometry, hydrodynamic method,
StateRedist controls, initial state, and final-time settings.

Generalize configured level-set construction to an arbitrary rectangular
region. Build root and fine geometry independently from the same plane or
circle definition, then require the established AMR measure compatibility
before PCM initialization.

At each root step, evaluate active-cell CFL limits on both synchronized levels
and choose `min(dt_c, r*dt_f)`, clipped to the remaining time. Run the qualified
two-level transaction, recompute stability, and stop only at the final time.
Write complete coarse and fine CSV files separately.

## Consequences

`pelef_reactive_eb_amr_2d` now executes a configured static hierarchy from input
through verified output. The two-level stability rule protects every fine
substep, and failures cannot increment time or step count.

The driver intentionally rejects chemistry and molecular transport until those
operators are composed transactionally on both levels. It does not tag or
regrid, support multiple patches or deeper levels, touch a physical domain
boundary, checkpoint, or distribute ownership with MPI.
