# Decision 0042: accept owner chemistry collectively

## Context

The `0.49.0` bridge assigns every patch one MPI owner, but no reacting-flow
operator uses that authority. Chemistry is the first suitable physics path to
distribute because each cell reactor is local: it needs no same-stage neighbor
flux, while still exercising real thermodynamics, kinetics, patch state,
hierarchy synchronization, and transactional failure handling.

## Decision

Before chemistry, synchronize every reactive patch from its owner and take an
identical backup on all ranks. Traverse levels and flattened patches in the
same order everywhere. Only the owner calls the existing serial chemistry
kernel. After each patch call, reduce the success flag with communicator-wide
logical AND. Broadcast an accepted owner's state, temperature, and ghost
storage before proceeding to the next patch.

After all patches are accepted, perform the existing deepest-to-root
average-down, recover temperatures, and rebuild physical, coarse/fine, and
sibling ghosts on every replica. If any owner or synchronization step fails,
restore the pre-chemistry backup on every rank and report one rejected
operation. Count only accepted owner calls, so a successful interval advances
each patch exactly once globally regardless of rank count.

## Consequences

MPI ownership now controls a real reacting-flow operator and preserves the
serial patch-tree result, conserved composite quantities, hierarchy
synchronization, and rollback semantics. The serial AMR module remains free of
MPI; only its chemistry and ghost-refresh operators become public building
blocks for the parallel scheduler.

This does not distribute flux-coupled operators. Hydro and molecular transport
still require owner-only recursive scheduling, stage/face communication,
shared time-integrated flux ownership, and distributed reflux. Replicated
storage also remains until migration and sparse local allocation are added.
