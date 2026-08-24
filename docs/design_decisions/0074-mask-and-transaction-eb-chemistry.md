# Decision 0074: mask and transact EB chemistry

## Context

The runnable EB application advances multispecies general-EOS hydro but rejects
chemistry. The existing 2D chemistry operator integrates every Cartesian cell
in place. Calling it directly would react covered cells and could leave earlier
cells mutated if a later reactor integration failed.

## Decision

Extend the shared 2D chemistry operator with an optional logical active mask.
When present, only true cells undergo conserved-to-primitive recovery and
constant-volume reactor integration. Perform all cell work on candidate state
and temperature arrays, committing them only after every selected cell
succeeds. Existing unmasked callers retain the same public behavior with a
stronger transaction guarantee.

Compose the EB step as reaction for half an interval, the complete hydro and
StateRedist transaction for one interval, then reaction for the second half.
Keep the caller's output equal to its input until the entire sequence succeeds.
Load the selected elementary or full H2/O2 mechanism in the EB executable just
as the regular reactive applications do.

## Consequences

Covered EB cells are bitwise inert under chemistry, active cells use the same
qualified reactor implementation as regular 2D, and any chemistry, Riemann,
redistribution, or EOS failure rolls back the full split step. Chemistry does
not add intercell communication, so this ownership contract can later be
applied independently on AMR patches. EB molecular transport remains separate.
