# 0224: Bind selected static EB AMR 2D restart to its generated context

## Status

Accepted for `0.232.0`.

## Context

The installed selected serial reactive EB AMR 2D application already passes
the configure-time bundle SHA-256, generated `explicit` or `implicit`
chemistry policy, and normalized bundle-order composition into the shared
static two-level driver. Its checkpoint schema 3 records the complete fixed
configuration, embedded-wall controls, patch, clock, diagnostics, and
coarse/fine state, but it cannot prove that a selected restart uses the same
generated mechanism context.

Changing schema 3 would break the frozen fixed checkpoint byte stream.
Accepting it from a selected reader, or treating partial selected metadata as
a fixed call, would silently remove the provenance gate.

## Decision

Keep two exclusive static two-level schemas in
`reactive_eb_amr_2d_driver_mod`:

- fixed calls omit selected metadata and retain schema 3 byte-for-byte;
- selected calls provide bundle SHA-256, integrator policy, and composition
  together and use schema 4;
- partial metadata, schema downgrade, and either cross-schema read fail.

Schema 4 writes `SELECTED_CONTEXT`, the canonical 64-digit lowercase bundle
SHA-256, policy, species count, and normalized bundle-order mole fractions
immediately after the common header. The existing species order,
configuration, embedded-wall identity, numerics, patch, clock, coarse/fine
state, and end marker follow unchanged. Context and all body records are read
into private candidates; public restart targets are assigned only after the
entire file validates.

The shared application selects fixed or selected persistence explicitly and
propagates actionable checkpoint failure context. The selected front end
permits scheduled writes, intentional stop, and restart only for static
exactly-two-level mode. Input, coarse output, fine output, checkpoint, and
restart aliases are rejected before the mechanism or any output is touched.

## Consequences

A selected static two-level checkpoint resumes only with the same generated
bundle, chemistry policy, normalized composition, species order, immutable
configuration, wall controls, patch, and field records. A separate-process
full-H2/O2 run reproduces uninterrupted coarse and fine CSV files byte-for-byte
and also matches the fixed uninterrupted reference.

The bundle SHA identifies provenance; it is not a digest of the checkpoint
payload. Finite physically consistent tampering, crash-atomic replacement,
schema conversion, selected dynamic/three-level/multipatch persistence,
MPI EB AMR restart, and scalable I/O remain separate work.
