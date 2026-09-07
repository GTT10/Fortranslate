# 0230: Dispatch selected mechanisms in sparse MPI EB patch-tree 2D

## Status

Accepted for `0.238.0`.

## Context

The fixed sparse MPI reactive EB patch-tree 2D application already owns a
rank-distributed, dynamically regridded hierarchy and composes chemistry,
transport, hydro, CSV output, and schema-8 restart. Configure-time selected
mechanisms already run in serial EB AMR 2D, but the sparse MPI EB chemistry
wrappers drop the generated `explicit` or `implicit` policy and there is no
selected MPI EB application target.

The fixed schema-8 checkpoint fingerprint does not identify a selected bundle,
bundle-order composition, or generated integrator. Reusing it for a selected
run would therefore allow restart under the wrong chemistry context. Adding a
selected persistent schema in the same increment would mix lifecycle dispatch
with a separate transactional-format change.

## Decision

Introduce one shared sparse MPI reactive EB patch-tree 2D application module.
The fixed and configure-time selected front ends load their own mechanism data
and call the same lifecycle. The selected call supplies a complete context:

- the 64-hex-digit configure-time bundle SHA-256;
- finite, nonnegative, normalized bundle-order mole fractions;
- the generated `explicit` or `implicit` chemistry integrator.

Every rank must present the same complete selected context before root-state
initialization. The selected composition initializes the root owner, and the
integrator is forwarded through both chemistry half-steps and the public
chemistry, full-physics, and to-time sparse MPI EB entry points. Invalid or
rank-dependent integrator policy is rejected before candidate publication.

The selected front end rejects input/output aliases and rejects every nonempty
`checkpoint_file` or `restart_file`, positive `checkpoint_interval`, or true
`checkpoint_stop_after_write` before loading the generated mechanism.
The shared application repeats the persistence rejection as a defensive
boundary. Fixed schema-8 checkpoint/restart behavior remains unchanged.

## Consequences

A pinned full-H2/O2 selected target can execute the owner-only root,
dynamically regrid the sparse EB patch tree, and produce rank-count-independent
composite output on one, two, and four ranks. Its one-rank result must match the
fixed full-H2/O2 reference exactly. Fixed checkpoint bytes remain covered by a
freeze gate.

This increment does not add selected MPI EB checkpoint/restart, three-level or
multipatch MPI execution, runtime mechanism loading, payload authentication,
crash-atomic or scalable I/O, CFD CVODE, MPI load balancing, performance or
thread qualification, detailed-fuel validation, or external PeleC field
parity. Selected MPI EB persistence is a separate schema decision.
