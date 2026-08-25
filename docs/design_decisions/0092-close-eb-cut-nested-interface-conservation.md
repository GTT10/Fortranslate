# Decision 0092: close EB-cut nested-interface conservation

## Context

The first recursive three-level hydro kernel rejected a finest coarse/fine
interface crossed by the embedded boundary. Reusing only the two-level flux
register on that geometry produced a repeatable mass and energy defect after
the second middle substep. The regular-interface path remained conservative,
showing that the recursive schedule and outer register were sound but the
nested cut-interface update required an additional multilevel closure.

## Decision

Record the authoritative middle/finest composite integral before every middle
interval. Integrate the signed Cartesian flux through the exterior boundary of
the middle mesh over that interval. After inner reflux and average-down,
compare the resulting composite integral with the pre-update integral plus
that boundary contribution.

Correct density, total energy, and every species. Do not correct momentum,
because a slip embedded boundary supplies a physical pressure force. Require
the species residual sum to close to the density residual within scaled
roundoff. Spread the conserved residual uniformly per fluid volume over active
middle cells outside the finest patch, then recover every recipient
temperature through the EOS. Publish no hierarchy if correction or recovery
fails.

## Consequences

The qualified plane-EB finest interface can cross cut and covered regions while
the three-level transaction retains composite mass, total energy, species,
positive thermodynamics, synchronization, and rollback. The recipient set is
outside the authoritative finest rectangle and therefore survives subsequent
deepest-first synchronization.

This is a global conservative multilevel closure. It does not claim parity
with PeleC's local multilevel EB redistribution or define an arbitrary-depth
recipient graph. Chemistry splitting, a public three-level lifecycle,
regridding, checkpointing, output, transport, and MPI remain separate work.
