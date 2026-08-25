# Decision 0086: connect the reactive EB patch-set application

## Context

The qualified reactive EB patch-set kernel could plan, rebuild, advance, and
roll back multiple separated fine rectangles, but the public executable still
constructed and owned the earlier single-patch hierarchy. Consequently users
could not select the patch-set lifecycle from an input file, exercise its
set-wide CFL and regrid operations, or write its child states. The formatted
checkpoint schema also has only one optional fine-patch payload.

## Decision

Add an explicit `multipatch_enabled` input switch and a maximum tag-gap control
to the reactive EB AMR configuration. Keep the existing single-patch path as
the default. When enabled, convert the configured seed rectangle into a valid
one-plan collection, build the patch set, and use collection planning for
initial and periodic temperature-driven regridding.

Select each coarse interval from the minimum root and refinement-scaled
all-child CFL limits. Advance it with the qualified set-wide Strang
chemistry/hydrodynamics transaction, and publish time, step, minimum-timestep,
and regrid counters only after the corresponding operation commits. Preserve
the existing empty-tag removal policy and treat an unchanged collection as a
topology no-op.

Write the synchronized root through the existing output path. Write each
active child in deterministic collection order by inserting
`_patchNNNN` before the configured fine CSV extension. Reject any checkpoint
or restart control in multipatch mode before hierarchy initialization rather
than silently serialize an incomplete collection.

## Consequences

The public serial EB AMR application can now follow two disconnected reactive
features with two conservative fine rectangles, including set-wide timestep
selection, regridding, physics advancement, diagnostics, and output. Existing
single-patch inputs and their checkpoint/restart behavior remain unchanged.

Multipatch checkpoints require a future schema revision. EB molecular
transport, deeper levels, physical-boundary fine patches, and distributed
patch ownership also remain outside this decision.
