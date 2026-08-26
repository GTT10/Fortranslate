# Decision 0150: Own the EB patch-tree stable-step clock

## Context

The arbitrary-depth numerical EB tree can select a hyperbolic CFL interval and
execute one atomic `R-T-H-T-R` step, but callers still have to repeat those
operations and manage time, accepted-step accounting, stop-time clipping, and
failure publication. Full physics also requires the explicit mixture
transport limit, which was not part of the hydro-only tree selector.

## Decision

Add one read-only selector that visits every active runtime node, evaluates
both the hyperbolic and enabled explicit transport limits, multiplies each
local interval by its cumulative ancestor refinement product, and publishes
the minimum root interval.

Compose that selector with the qualified full-physics transaction in a public
target-time loop. Clip each interval to the remaining time, execute the step on
a private tree candidate, and commit the tree, time, step count, minimum
accepted interval, minimum limiter theta, and per-level physics counts only
after success. Treat the caller's maximum step count as a termination bound,
not permission to expose a failed candidate.

## Consequences

Callers can now advance a static arbitrary-depth or branching EB numerical
tree to an exact target time without reproducing stability or accounting
logic. The accepted state always matches the published clock. A first-step
failure is an exact rollback; a later failure or exhausted step budget retains
only the successful prefix.

This decision does not schedule topology rebuilds, write or restore
checkpoints, or distribute arbitrary-depth tree fields across MPI ranks.
