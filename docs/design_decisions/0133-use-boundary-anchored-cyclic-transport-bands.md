# Decision 0133: use boundary-anchored cyclic transport bands

## Context

Sparse transport root Euler stages use a six-row dependency guard around each
owned target tile. In `0.140.0`, any target touching a periodic y boundary used
the complete root because an ordinary clipped band would make the serial
kernel wrap its two ends onto unrelated interior rows.

A cyclic band can be assembled from rows near both global y boundaries, but a
naive ordering introduces a second, artificial adjacency where the two selected
global ranges meet. Transport and StateRedist run over the whole temporary
geometry, so a cell on the required guard must not depend on that seam.

## Decision

Represent a periodic-edge band as two increasing, globally contiguous source-
row fragments. Keep global row one at the lower temporary boundary and global
row `ny` at the upper temporary boundary, preserving the physical periodic wrap
at the serial kernel's outer ends. Place the target at the corresponding outer
edge and retain its qualified six-row dependency guard.

Add one seam-isolation row beyond the six-row guard on both fragments. Use the
cyclic form only when the target plus both protected fragments is smaller than
the root; otherwise retain the complete-root fallback. Split assembly at both
source-owner boundaries and the deliberate global-row gap, and count every
point-to-point fragment and computed band cell exactly.

Copy EB metric arrays by selected global row. Shift absolute cut-boundary
centroid y coordinates by the compact row displacement, while leaving
non-cut sentinel values and local centroid/normal metrics unchanged. Preserve
the physical lower and upper y-face metrics at the band ends.

## Consequences

Periodic edge tiles on sufficiently tall roots no longer require an
unconditional complete-root input band. The dedicated 14-by-21 gate exercises
finite cyclic bands at four and eight ranks and retains complete-root fallback
at one and two ranks.

The target loop remains deterministic and sequential, so no wall-time speedup
is claimed. A complete post-compute root bundle also remains for fine-child
exterior construction, flux registers, and reflux. Those interfaces are a
separate distributed-data boundary.
