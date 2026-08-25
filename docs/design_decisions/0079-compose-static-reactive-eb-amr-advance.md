# Decision 0079: compose a static reactive EB AMR advance

## Context

The static EB hierarchy can already validate one aligned fine rectangle,
average down by fluid volume, accumulate open-area flux mismatch, and re-reflux
cut interface cells. Those kernels did not initialize the fine state, supply
coarse data at a fine-patch boundary, or schedule the level advances over one
shared physical interval.

## Decision

Add a serial reactive hydrodynamic composition for one strictly internal fine
rectangle. Initialize it by piecewise-constant parent injection and recover
temperature on each active child. Advance the coarse level once over `dt`, then
advance the fine level `r` times over `dt/r`.

Represent the fine-patch exterior with four face-indexed state and temperature
arrays. On every open boundary face, interpolate the adjacent coarse conserved
state between the beginning and accepted end of its step, then recover
temperature with the reactive EOS. Use the start of a fine interval for PCM
and its midpoint for the time-centered characteristic-PLM trace. Keep the
coarse-to-fine spatial transfer piecewise constant for this milestone.

Accumulate the exact centroid fluxes used by each level into the EB register.
After all fine substeps, apply reactive re-reflux and then reactive EB
average-down. Hold all caller outputs at the original hierarchy until the
complete chain succeeds.

## Consequences

One coarse and `r` fine hyperbolic steps now share consistent coarse-time
boundary data and conservative synchronization. Uniform stationary reactive
states survive the diagonal EB hierarchy to roundoff, including covered cells
and the composite conserved integral. Any invalid solver, boundary state, EOS
conversion, reflux, or restriction rejects the complete interval.

The patch must remain static, rectangular, aligned, and strictly inside the
coarse domain. This does not add spatially sloped coarse interpolation,
physical-domain-touching patches, chemistry or transport composition, dynamic
regridding, multiple patches, deeper levels, or MPI ownership.
