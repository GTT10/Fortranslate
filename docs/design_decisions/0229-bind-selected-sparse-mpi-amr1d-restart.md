# 0229: Bind selected sparse MPI AMR 1D restart to its context

## Status

Accepted for `0.237.0`.

## Context

The fixed sparse MPI reactive AMR 1D application already distributes a
dynamic patch tree and can rebuild ownership for the current communicator.
Its checkpoint uses
`PELEF_AMR_PATCH_TREE_REACTIVE_1D_CHECKPOINT` schema 1 and deliberately stores
neither an MPI rank count nor an owner map. The selected-mechanism application
did not exist, the sparse chemistry wrappers dropped the generated integrator
policy, and schema 1 could not identify the selected bundle, composition, or
original conservation baseline.

Reusing schema 1 for selected runs would permit restart under the wrong
generated mechanism. Persisting an owner map would make restart depend on the
writer's communicator rather than on the physical hierarchy. Modifying schema
1 would also invalidate established fixed checkpoints.

## Decision

Keep fixed magic/schema 1 byte-for-byte. A selected sparse MPI AMR call uses
schema 2 under the same magic. Immediately after the species order it writes:

- `SELECTED_CONTEXT` with the configure-time bundle SHA-256, generated
  `explicit` or `implicit` integrator, species count, and finite normalized
  bundle-order mole fractions;
- `COMPOSITE_BASELINE` with the finite initial composite conserved integrals.

The remaining body stores the rank-neutral patch tree, committed fields,
clock, regrid diagnostics, and terminal marker. The selected reader validates
the expected context, baseline closure, complete payload, recovered physical
state, terminal marker, and end-of-stream in private candidates before
publishing any caller target. Failed reads and invalid writes are
non-destructive. Fixed and selected readers reject the other schema.

One shared MPI application module owns both fixed and selected lifecycles.
Each selected rank must present the same complete context before state
initialization or output. Bundle-order composition initializes the root and
the generated integrator is forwarded through every sparse chemistry half-
step. Restart restores the original baseline, rebuilds distribution for the
current communicator, and never consumes persisted ownership metadata. The
selected front end rejects lexical aliases among input, output, checkpoint,
and restart paths before loading the generated mechanism.

## Consequences

A pinned full-H2/O2 run stops after its first coarse step on one rank and
continues independently on two and four ranks. Both continuations are
byte-identical to uninterrupted one-, two-, and four-rank selected runs and to
the fixed one-rank reference. Rank-dependent context, bundle, integrator,
composition, species-order, and baseline mutations are rejected without an
output artifact. An invalid integrator fails collectively at the first sparse
chemistry half-step with the solution unchanged.

This is selected configure-time dispatch and rank-neutral formatted restart
for the one-dimensional sparse MPI AMR patch tree. It does not add selected
MPI EB AMR, runtime mechanism loading, crash-atomic or scalable I/O, payload
authentication, thread safety, performance qualification, detailed-fuel
validation, or external PeleC field parity.
