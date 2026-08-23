# Decision 0019: transactional MPI reactive 1D composition

## Status

Accepted for PeleF `0.24.0`.

## Context

The serial full-H2/O2 path already separates cell-local implicit chemistry,
general-EOS molecular transport, and conservative hydro. Distributed execution
must preserve that operator ordering while preventing one rank from accepting a
trial that another rank rejects.

## Decision

Use uneven contiguous one-dimensional blocks with one periodic ghost cell on
each side. Keep conserved state and derived temperature rank-local, exchange
both fields before face operators, and use communicator-wide reductions for
timestep limits, success flags, conservation diagnostics, and output assembly.

Compose a distributed step as

```text
chemistry(dt/2)
transport(dt/2)
hydro(dt)
transport(dt/2)
chemistry(dt/2)
```

on a trial copy. If any rank reports failure, reduce that decision across the
communicator, restore the complete pre-trial state on every rank, and retry at
half the interval. A successful attempt is committed only after all ranks agree.

## Consequences

- rank-count changes follow the same accept/reject history;
- local state is not replicated merely to simplify chemistry or transport;
- the existing serial physics kernels remain the reference implementation;
- ordered root gather is confined to verification output, not timestepping;
- the current claim is one-dimensional uniform-grid MPI, not AMR or
  multidimensional load balancing.

## Verification

Debug and Release CI run the halo, multispecies, implicit chemistry, molecular
transport, and coupled reactive gates with 1, 2, 4, and 8 ranks. Complete
gathered fields must agree within `5e-13` relative tolerance, and the
MPI-enabled build must retain the complete serial regression suite.
