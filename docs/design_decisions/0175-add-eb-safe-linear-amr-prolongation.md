# Decision 0175: Add EB-safe limited-linear AMR prolongation

## Context

Reactive 2D EB regrids initialize every new fine child with its parent state.
That PCM contract is conservative and robust across cut geometry, but it
discards smooth coarse-level variation even when a complete parent/child block
is regular. Public orchestration should not change until a higher-order kernel
has an isolated conservation and admissibility gate.

## Decision

Add a separate low-level prolongation routine using component-wise MC slopes
of the conserved state. Enable those slopes only for a regular coarse parent
whose complete fine-child block is regular. The Cartesian child offsets have
zero mean, preserving the parent state under average-down. Recover each child
temperature through the EOS instead of interpolating temperature.

Keep PCM for cut, covered, or topology-mismatched parents. If any reconstructed
regular child is not EOS-admissible, discard every linear child of that parent
and retry the parent with PCM. Publish neither state nor temperature until the
complete patch succeeds.

## Consequences

Library callers gain smooth second-order initialization away from the EB while
the qualified cut-cell behavior remains unchanged. The new kernel is exactly
conservative on accepted Cartesian parents and transactionally preserves the
existing admissibility boundary. Public namelist selection and checkpoint
identity remain separate work, so the current application default stays PCM.
