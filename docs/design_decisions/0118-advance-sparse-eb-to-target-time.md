# Decision 0118: advance sparse EB owners to a target time

## Context

The sparse MPI reactive EB hierarchy can select one owner-local stable coarse
interval and execute one direct `R-T-H-T-R` transaction, but callers still
have to own the public clock, repeat timestep selection, clip the final
interval, and reconcile diagnostic publication with transactional failures.

## Decision

Add a public collective loop that accepts an existing time and total step
count, a target time, and a total-step limit. Require communicator consensus
for those controls and for every numerical tolerance and physics switch. Before
each step, select a fresh sparse hydro/transport stability limit, clip it to the
remaining target interval, and execute the existing direct sparse
`R-T-H-T-R` transaction.

Publish time, total steps, minimum accepted dt, owner operator counts, limiter
minimum, and timestep root-gather traffic only after the complete split step
commits. If a later selection or physics stage fails, or the total-step limit
is reached, retain every earlier committed state and publish exactly the
diagnostics belonging to that prefix. A request already at its target is a
successful zero-step operation.

## Consequences

Static sparse MPI EB AMR now owns a qualified full-physics clock without
materializing fine fields. The API has explicit restart-friendly time and step
inputs and exposes committed work separately from the cumulative total step.

This loop intentionally keeps topology fixed. Dynamic sparse regridding,
checkpoint/restart, parallel output, and decomposition of the level-wide root
physics kernel remain future lifecycle work.
