# Decision 0083: compose chemistry on reactive EB AMR

## Context

The serial EB AMR application had a qualified lifecycle for one fine rectangle,
but it rejected chemistry even though the single-level EB driver and the
regular 1D AMR driver already had transactional Strang operators. Running a
reactive EB hierarchy therefore stopped at hydrodynamics.

## Decision

Wrap the existing two-level EB hydro/subcycling/reflux transaction with an
active-cell reaction half-step on each level. After the second reaction
half-step, average the fine reactive state into its covered coarse parents so
the public hierarchy is synchronized at the end of the interval.

Keep all intermediate coarse and fine state and temperature arrays private.
Publish them only when both chemistry halves, hydrodynamics, EOS recovery,
reflux, and average-down succeed. Build a separate active mask from each EB
geometry so covered cells never enter chemistry. When the lifecycle has no
fine patch, reuse the existing single-level EB Strang operator.

Load the elementary or full H2/O2 reaction mechanism in the public EB AMR
application through the same input-driven path as the regular and single-level
EB applications. Continue rejecting molecular transport before any timestep is
accepted.

## Consequences

One serial, ratio-aligned, strictly internal EB fine rectangle can now react
while active, collapse conservatively, and continue reacting on the root. A
failure after chemistry begins cannot publish a partially reacted hierarchy.

Molecular transport, multiple simultaneous patches, deeper levels,
checkpoint/restart, and distributed EB ownership remain outside this decision.
