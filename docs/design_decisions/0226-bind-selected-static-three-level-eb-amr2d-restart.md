# 0226: Bind selected static three-level EB AMR 2D restart to its context

## Status

Accepted for `0.234.0`.

## Context

The fixed serial three-level reactive EB AMR 2D lifecycle already has two
released checkpoint contracts. Static calls use the
`PELEF_REACTIVE_EB_AMR_THREE_LEVEL_2D_CHECKPOINT` magic with schema 3, while
dynamic-finest calls use a separate magic with schema 4. Both byte streams are
frozen. Selected persistence was qualified only for the two-level topology:
schema 4 for static calls and schema 5 for dynamic calls.

Reusing the fixed static schema 3 would omit the generated mechanism identity.
Reusing the fixed dynamic magic would misstate the topology class. A selected
three-level restart also needs the original composite-integral baseline so a
resumed cumulative conservation diagnostic has the same reference as an
uninterrupted run.

## Decision

Keep both fixed three-level byte streams unchanged. A selected static
three-level call uses the existing static magic with schema 4. Immediately
after the header it writes:

- an all-or-none `SELECTED_CONTEXT` block containing the configure-time bundle
  SHA-256, generated `explicit` or `implicit` integrator, species count, and
  finite normalized bundle-order composition;
- a `COMPOSITE_BASELINE` block containing the finite `nvar`-component initial
  three-level composite integrals.

The existing static three-level body then records species, geometry, wall and
numerical policies, both nested patches, clock and diagnostics, root/middle/
finest state and temperature fields, and the terminal marker. The reader
derives the expected schema from the fixed/selected call contract. It parses
context, baseline, topology, every field, recovered-temperature consistency,
and the terminal marker into private candidates before publishing any caller
target.

The selected application admits this envelope only for a serial static
three-level hierarchy with exactly one patch per level. Selected dynamic
three-level, multipatch, and dynamic-parent modes remain rejected. Lexical
alias checks include input, root, middle, finest, checkpoint, and restart
paths before a file is touched.

## Consequences

A pinned full-H2/O2 process can stop after its first coarse step and resume in
another process to byte-exact uninterrupted root, middle, and finest CSV
files. The selected uninterrupted files also equal the fixed static reference
at every level. Restoring the persisted baseline makes the final cumulative
conservation diagnostic textually identical. Changed context, corrupt
baseline, fixed/selected cross-schema input, truncation, or any invalid field
publishes no candidate; every output target remains in its documented
empty/default failure state.

The selected schema records provenance but does not authenticate its payload.
Selected dynamic three-level, multipatch, dynamic-parent, MPI EB AMR,
crash-atomic or scalable I/O, schema migration, and payload authentication
remain separate work.
