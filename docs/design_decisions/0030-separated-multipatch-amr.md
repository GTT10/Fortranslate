# Decision 0030: qualify separated patches before same-level box exchange

## Context

The original AMR data model stores one contiguous fine array per level.
Representing several boxes introduces two distinct interfaces: coarse/fine
sides around separated patches and fine/fine sides where boxes are adjacent.
Treating an adjacent fine/fine side as two coarse/fine sides would apply
incorrect reflux corrections and duplicate its interface flux.

## Decision

The first multipatch type stores ordered parent-index intervals and requires at
least one uncovered parent cell between consecutive patches. Tag clustering
coalesces candidate bounds that touch or overlap before constructing the set.
Each remaining patch therefore owns zero, one, or two genuine coarse/fine
sides and can reuse the qualified prolongation, ghost-fill, flux-register,
reflux, and average-down kernels.

Set-wide operations are transactional. Composite integration starts from the
complete parent integral, removes every disjoint covered interval, and adds
each fine patch. Regridding synchronizes all old patches first and copies every
same-resolution old/new fine intersection after conservative prolongation.

## Consequences

Fixed two-level hydro can advance several separated patches conservatively
without inventing same-level communication. Adjacent boxes must be coalesced,
which can refine extra cells. Supporting independently owned adjacent boxes
later requires explicit same-level ghost exchange, interface ownership, and a
single flux per fine/fine face. At the `0.38.0` milestone chemistry, transport,
arbitrary-depth patch trees, and MPI distribution remained separate;
Decision 0031 subsequently integrates fixed two-level chemistry and transport.
