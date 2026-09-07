# 0222: Dispatch selected mechanisms in serial EB AMR 2D flow

## Status

Accepted for `0.229.0`.

## Context

The serial reactive EB AMR 2D application already owns a static two-level
hierarchy, fine subcycling, EB reflux, conservative average-down, chemistry,
molecular transport, configured embedded walls, deterministic level output,
and diagnostics. Its program and driver nevertheless supplied the committed
full-H2/O2 composition and chemistry policy directly. Linking a generated
mechanism to `pelef_core` would also retain the fixed mechanism namespaces,
while copying the AMR loop would create a second hierarchy algorithm.

The same driver also exposes dynamic, three-level, multipatch, and restart
paths. The existing checkpoint schema does not record the selected bundle
fingerprint or generated integrator policy, so those modes cannot be claimed
by merely accepting dynamic species arrays at startup.

## Decision

Extract `reactive_eb_amr_2d_application_mod`. Fixed and selected front ends
load their own thermo, reaction, transport, and composition data and then call
one application driver for hierarchy construction, subcycled advancement,
output, invariants, and diagnostics. Thread an optional bundle-order
composition through regular and configured boundary initialization. Thread an
optional generated chemistry-integrator policy through both coarse and fine
chemistry half-steps. An absent argument preserves the fixed front end.

Build `pelef_selected_reactive_eb_amr_2d_runtime` in a private Fortran module
directory above the selected single-level EB runtime. Configure and install
`pelef_reactive_eb_amr_2d_selected` from each complete schema-1 bundle without
linking `pelef_core` or either committed mechanism.

Qualify only the static two-level serial path. Reject dynamic regridding,
three-level and multipatch operation, checkpoint creation, and restart before
model loading or output creation. Also require input, coarse-output, and
fine-output paths to be lexically distinct. Validate exact-name composition
and every requested initial, hotspot, isothermal-wall, and reflected-wall
temperature against the selected mechanism's common NASA7 interval.

## Consequences

Selected two- through 32-species bundles can execute the established static
two-level 2D EB hierarchy with native explicit or implicit chemistry and
molecular transport. A shuffled two-species fixture requires positive,
closed, deterministic coarse and fine output. A shuffled pinned full-H2/O2
case requires active chemistry on both levels and byte identity with the fixed
front end.

This milestone does not qualify dynamic, three-level, or multipatch selected
EB AMR, checkpoint/restart, selected EB 3D, MPI EB AMR, runtime mechanism
loading, CVODE-backed CFD chemistry, mechanisms above 32 species, thread
safety, performance, detailed fuels, external PeleC field parity, or
experiment validation.
