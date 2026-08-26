# Decision 0141: Keep root transport results on tile owners

## Context

Direct state, coarse-flux, and corrected-support routing removed every child
dependency on the temporary complete root transport result. Sparse transport
still sent each remote tile's stage-start state, stage-end state, temperature,
and flux rows to the root physics owner. Only physical-boundary flux closure
continued to read that assembled bundle.

## Decision

Keep transport Euler results exclusively in owner-local tile state and flux
records. Remove the post-compute tile-result send and every complete root
transport result allocation. Compute left and right physical-boundary flux
contributions on every tile owner, lower and upper contributions on the edge
tile owners, and combine the resulting conserved `nvar` vector with
`MPI_Allreduce` before the established distributed closure.

Count only point-to-point halo and direct child-fragment messages as sparse
root transport traffic. Preserve ordered child correction, final tile-local
publication, rollback, and the explicit materialization boundaries used by
checkpoint/output compatibility. Hydro is unchanged.

## Consequences

Sparse molecular transport has no post-compute complete-root numerical field
or root-physics-owner result bottleneck. Cut-interface closure adds one compact
communicator sum per Euler stage that actually contains a cut interface.
Qualification must cover serial field and conservation parity, exact traffic,
root-only cyclic bands, owner work, and failure rollback in Debug and Release
at one, two, four, and eight ranks. This decision does not claim a measured
speedup.
