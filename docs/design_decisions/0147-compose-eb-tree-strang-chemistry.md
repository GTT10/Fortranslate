# Decision 0147: Compose EB patch-tree Strang chemistry

## Context

The `0.154.0` numerical EB tree advances hydrodynamics at arbitrary runtime
depth, but active-cell chemistry composition remains limited to fixed-depth and
sibling patch-set APIs. Invoking those APIs per relation would repeat ancestor
reaction updates and would not provide one rollback boundary for a branching
tree.

## Decision

Apply the existing 2D chemistry integrator once to every runtime patch with an
active mask derived from that patch's EB geometry. A standalone public
chemistry transaction synchronizes the candidate deepest first and publishes
optional per-level patch-call counts only after every patch succeeds.

Compose the full reaction/hydrodynamics interval on one private tree candidate:
apply the first reaction half-step to every patch, invoke the qualified
recursive tree hydrodynamics transaction, apply the second reaction half-step,
then synchronize deepest first. Do not publish state, temperature, chemistry
counts, or hydro counts until the final candidate validates. Chemistry-disabled
use retains the same entrypoint and reduces to the qualified hydro operation.

## Consequences

Arbitrary-depth and branching EB trees now support active-cell reaction and
transactional `R-H-R` Strang splitting without a fixed level count. A
four-level two-branch gate verifies exact chemistry and hydro schedules,
conservation, thermodynamics, and rollback after chemistry has already changed
the private candidate. A runtime three-level chain retains field and
temperature parity with the established fixed-depth Strang path.

This decision does not add molecular transport, a public time loop, dynamic
tagging, checkpoint I/O, or MPI ownership for the numerical tree.
