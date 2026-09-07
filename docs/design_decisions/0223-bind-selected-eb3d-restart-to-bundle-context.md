# 0223: Bind selected EB 3D restart to its generated context

## Status

Accepted for `0.231.0`.

## Context

The selected serial planar reactive EB 3D front end already supplies a
configure-time bundle SHA-256, a generated `explicit` or `implicit` chemistry
policy, and a normalized bundle-order composition to the shared application.
The existing checkpoint stores complete thermo, reaction, transport,
configuration, geometry, state, temperature, and cumulative diagnostics, but
schema 1 has no selected-context record. Accepting it from a selected front
end could therefore continue the same numerical records under a different
generated policy or initial-composition contract.

Changing every checkpoint to schema 2 would also break the frozen fixed
schema-1 byte stream. Falling back from a selected reader to schema 1, or from
an incomplete selected writer to schema 1, would silently remove the identity
gate that this increment is intended to add.

## Decision

Keep two exclusive schemas in `reactive_eb_3d_checkpoint_mod`:

- a call with no selected metadata writes and reads fixed schema 1 exactly as
  before;
- a call with bundle SHA-256, integrator policy, and composition writes and
  reads selected schema 2;
- partial metadata, schema downgrade, and cross-runtime schema reads fail.

Schema 2 places `SELECTED_CONTEXT`, the canonical 64-digit lowercase bundle
SHA-256, the generated integrator, species count, and normalized
bundle-order mole fractions immediately after the header. The reader validates
all context fields before parsing private geometry/state candidates and keeps
the existing publish-only-after-end-marker transaction.

The shared application explicitly selects the fixed or selected checkpoint
API. The selected front end enables scheduled checkpoints, intentional stop,
and restart only after the complete context is available. It also rejects
lexical input/checkpoint and input/restart aliases before any file can replace
the input.

## Consequences

A selected checkpoint can resume only with the same generated bundle,
chemistry policy, bundle-order composition, complete model records, geometry,
and immutable numerical controls. A two-process selected elementary run
reproduces the uninterrupted final CSV byte-for-byte and retains cumulative
diagnostics. Fixed schema-1 checkpoint and stopped-output hashes remain
unchanged.

The stored bundle SHA identifies the mechanism bundle; it is not a digest of
the checkpoint payload. Truncation, malformed records, nonfinite or physically
inconsistent state, and context mismatch are rejected, but a finite,
physically consistent payload alteration is not cryptographically detected.
Crash-atomic replacement, schema conversion, scalable or distributed I/O,
general 3D EB geometry, AMR/MPI EB 3D restart, and tamper evidence remain
separate work.
