# 0227: Bind selected dynamic three-level EB AMR 2D restart to its context

## Status

Accepted for `0.235.0`.

## Context

The fixed serial dynamic three-level reactive EB AMR 2D lifecycle already
stores dynamic-finest and dynamic-parent policy, the committed middle and
finest patches, regrid history, and all three fields under
`PELEF_REACTIVE_EB_AMR_DYNAMIC_THREE_LEVEL_2D_CHECKPOINT` schema 4. That byte
stream is frozen. Selected static three-level persistence uses the distinct
static magic with schema 4, but selected dynamic three-level calls were still
rejected.

Reusing fixed schema 4 would omit the generated mechanism identity. Reusing
the static selected format would misstate the topology class. A selected
restart also needs the original composite-integral baseline so resumed
cumulative conservation diagnostics retain the uninterrupted reference.

## Decision

Keep the fixed dynamic magic/schema 4 byte-for-byte. A selected dynamic
three-level call uses schema 5 under the same dynamic magic. Immediately after
the header it writes:

- an all-or-none `SELECTED_CONTEXT` block containing the configure-time bundle
  SHA-256, generated `explicit` or `implicit` integrator, species count, and
  finite normalized bundle-order composition;
- a `COMPOSITE_BASELINE` block containing the finite `nvar`-component initial
  three-level composite integrals.

The existing dynamic schema-4 body then records species, geometry, wall and
numerical policies, dynamic controls, committed middle and finest patches,
clock and regrid diagnostics, root/middle/finest state and temperature fields,
and the terminal marker. The reader derives the expected schema from the
fixed/selected call contract. It parses every record and checks recovered
temperatures and end-of-stream into private candidates before publishing any
caller target.

The selected application admits this envelope only for a serial exactly-
three-level hierarchy with one patch per level. Both the established finest-
only and dynamic-parent lifecycles are allowed. Selected multipatch mode
remains rejected. Lexical alias checks cover input, root, middle, finest,
checkpoint, and restart paths before a file is touched.

## Consequences

A pinned full-H2/O2 process moves both nested patches, stops after its first
coarse step, and resumes in another process to byte-exact uninterrupted root,
middle, and finest CSV files. The selected uninterrupted files also equal the
fixed dynamic reference at every level. Restoring the persisted baseline makes
the final cumulative conservation diagnostic textually identical. The
existing fixed dynamic regression independently pins the schema-4 checkpoint
SHA-256 so this increment cannot silently rewrite the frozen byte stream.

Changed context, corrupt baseline, fixed/selected cross-schema input,
static/dynamic cross-magic input, corrupt controls, patch, field, or terminal
record, truncation, or trailing content publishes no candidate. An invalid
write leaves an existing checkpoint untouched.

The selected schema records provenance but does not authenticate its payload.
Selected multipatch and MPI EB AMR persistence, crash-atomic or scalable I/O,
schema migration, and payload authentication remain separate work.
