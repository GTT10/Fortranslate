# Decision 0165: Expose the EB tree as a dedicated application

## Context

The arbitrary-depth 2D EB numerical tree had qualified initialization,
physics, dynamic topology, checkpoint/restart, and composite output APIs, but
only unit and MPI verification programs composed them. The established public
EB AMR application already branches among single-patch, sibling-multipatch,
and fixed three-level storage contracts.

## Decision

Add a separate installed `pelef_reactive_eb_patch_tree_2d` executable. Reuse
the existing `&reactive_2d`, `&embedded_boundary`, and `&eb_amr` namelists and
add one bounded `patch_tree_maximum_levels` control. Let a dedicated driver
initialize or restart the tree, apply initial and scheduled recursive tags,
select and commit full-physics root steps, invoke checkpoint cadence, and write
the qualified composite CSV.

Construct child EB geometry through an internal callback over the configured
plane or circle. Keep lifecycle counters in the driver and numerical topology
and fields in the patch-tree core.

## Consequences

Users can now run the arbitrary-depth serial EB path from an input file without
changing the behavior or output contracts of any fixed-depth mode. The new
executable can evolve independently toward application checkpoint/restart
parity and a sparse MPI counterpart. Input fields irrelevant to patch trees
remain accepted for compatibility but are not interpreted as static children.
