# Decision 0100: add three-level EB AMR transport recursion

## Context

The two-level EB AMR transport transaction already provides fine subcycling,
time-interpolated exterior states, diffusive reflux, reactive average-down,
and SSPRK2 rollback. A three-level hierarchy needs the same contract at two
nested interfaces while retaining a single authoritative state for each
physical fluid volume.

## Decision

Advance the root transport stage once. Advance the middle level once per root
substep and, within every middle interval, advance the finest level once per
middle substep. Give each parent/child interface an independent diffusive flux
register. Complete finest reflux, EB-cut conservation closure, and
middle/finest average-down before accumulating the middle result into the root
transaction. Finish with outer reflux, deepest-first three-level average-down,
and the same EB-cut closure at the outer interface when required.

Compose two complete synchronized Euler transactions as SSPRK2. Convert the
middle and finest parabolic limits to root time by the cumulative refinement
ratio, and insert the operator into both transport half-steps of the existing
three-level reactive Strang composition.

## Consequences

The public static and movable-finest three-level lifecycles can now run the
qualified molecular-transport subset while preserving composite conservation
and all-level rollback. The transport checkpoint restriction remains because
the existing schemas do not record transport compatibility.

This milestone does not add sibling-patch transport, coarse-to-fine spatial
slopes, thermal or catalytic embedded walls, distributed ownership, parallel
flux registers, or transport checkpoint/restart.
