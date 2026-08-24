# Decision 0046: establish rank-local MPI AMR patch storage

## Context

The `0.53.0` MPI AMR path executes physics only on patch owners, but every rank
still allocates every patch payload. A direct sparse physics implementation
needs a storage contract before halo, flux-register, and regrid communication
can be changed safely.

## Decision

Add a sparse reactive patch-tree type that replicates hierarchy, ownership,
and bookkeeping metadata while allocating state, temperature, and four ghost
arrays only for locally owned patches. Require exact agreement between local
allocation flags and the owner map, and require local patch and cell counts to
match distribution metadata.

Provide explicit transition operations. Scatter resolves an
owner-authoritative replica before retaining local payloads. Gather restores
local owner data into a supplied replica and synchronizes every field and
global counter. A same-hierarchy migration sends each complete payload from
its old owner and retains it only on its new owner.

## Consequences

Persistent reactive field storage is no longer multiplied by MPI rank count,
and ownership changes can preserve the exact solution without returning to a
persistent replicated representation. The transition migration currently
uses one patch broadcast at a time, so it is a correctness boundary rather
than the final scalable schedule.

The `0.53.0` physics operators are unchanged and still require replicas.
Direct physics on sparse storage, point-to-point halos and fluxes, and
topology-changing regrid migration remain separate work.
