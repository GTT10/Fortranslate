# Decision 0129: decompose sparse EB root hydro with targeted halos

## Context

The sparse MPI EB hierarchy stores each root row tile only on its owner, but
direct hydro gathers every tile to the first tile's owner before advancing the
complete root. That avoids all-rank replication yet leaves both root compute
and complete input storage concentrated on one rank. The replicated path has
already qualified a six-row owner-tiled decomposition against the serial EB
hydro transaction.

## Decision

For every root target tile, identify the source tiles intersecting its owned
rows plus at most six guard rows. A remote source owner sends one packed state
and temperature fragment directly to the target owner; same-owner fragments
are copied locally. The target extracts the matching EB geometry band and
invokes the established reactive EB level kernel.

The target publishes only its owned input rows, updated state and temperature,
x-face rows, and uniquely assigned y-faces. Remote target owners send one
packed result to the root owner. That owner assembles the temporary complete
start, result, and flux arrays expected by existing fine-child exterior,
reflux, and final sparse-scatter code. Count only actual point-to-point sends,
and defer public transfer, tile-advance, and computed-cell counters until the
whole sparse hydro transaction commits.

## Consequences

No complete sparse root input is assembled before hydrodynamic work, and each
root tile is advanced on its existing owner with bounded overlapping work.
Only owners of intersecting halo rows communicate; unrelated ranks allocate no
complete root fields. Serial numerical behavior and the established child
transaction remain unchanged.

The root owner still assembles complete temporary arrays after root compute to
serve level-wide child exterior interpolation and deterministic reflux. Sparse
transport, timestep selection, distributed flux-register interfaces, and
elimination of that post-compute compatibility bundle remain later work.
