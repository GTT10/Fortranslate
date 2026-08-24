# Decision 0064: materialize sparse AMR only at run boundaries

## Context

The owner-local sparse kernels now provide tagging, topology rebuilds, stable
timesteps, chemistry, hydro, and molecular transport, but users could not run
them from a namelist. The existing gather API also required a preallocated
replicated hierarchy, which is not available after dynamic sparse regridding.

## Decision

Add a public MPI application that initializes the configured root solution,
scatters it to deterministic owners, creates the initial tagged hierarchy on
those owners, and then performs stop-time-clipped sparse `R-T-H-T-R` advances.
At each configured cadence it invokes the owner-local tagged regrid and adopts
the rebuilt ownership map.

Expose `amr_mpi_work_exponent` in the reactive namelist. Values zero, one, and
two select cell, hyperbolic-subcycle, and parabolic-subcycle work estimates.
Patch-tree tagging no longer depends on the legacy two-level multipatch flag.

Add a materialization helper that allocates a replicated tree from sparse
hierarchy metadata and fills it from the authoritative owners. Use it only
after the time loop for composite integrals and ordered CSV output. The output
writer recursively replaces every covered parent interval with its children.

## Consequences

Normal evolution and topology changes retain globally single-copy fields while
the public program supports arbitrary configured AMR depth. Initial root setup
and final output use replicated storage, so this is not yet a scalable
checkpoint or parallel-I/O design. Restart, asynchronous output, and measured
runtime load balancing remain separate capabilities.
