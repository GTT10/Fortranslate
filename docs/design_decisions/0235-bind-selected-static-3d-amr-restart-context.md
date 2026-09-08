# 0235: Bind selected static 3D AMR restart to its complete context

## Status

Accepted for `0.243.0`.

## Context

The `0.242.0` serial and rank-local sparse-MPI static 3D AMR applications
dispatch configure-time selected mechanisms, but reject checkpoint and restart.
The retained fixed schema 2 contains both levels and a NASA7/numerical
fingerprint, but it does not identify the selected bundle, generated chemistry
integrator, bundle-order composition, ordered reaction mechanism, or chemistry
tolerances. Reading that schema for selected chemistry could therefore publish
state under a different source operator.

## Decision

Keep fixed checkpoint calls byte-identical at schema 2 and give selected calls
an exclusive schema 3 under the established static-3D-AMR magic.

- Selected mode requires all four context arguments: bundle SHA-256, generated
  integrator, normalized bundle-order mole fractions, and ordered reactions.
  Partial context is rejected before the target file is opened.
- Schema 3 records the SHA and integrator, composition, and every ordered
  reaction equation, kind, reversibility flag, reactant/product stoichiometry,
  forward/low/high Arrhenius records, Troe data, and third-body efficiencies.
  It also records the chemistry relative and absolute tolerances.
- The established body then binds ordered NASA7 records, thermodynamic model,
  solver, boundary condition, reconstruction, limiter, mesh and patch layout,
  CFL, chemistry/transport flags, clock, step, original composite baseline,
  reflux history, both conserved levels, recovered temperatures, terminal
  marker, and strict end-of-stream.
- The selected reader validates context and the complete body into private
  candidates. It publishes no caller state, time, counters, baseline, or
  diagnostics unless all checks and the final close succeed.
- MPI performs exact selected-context consensus before root I/O. Root reads or
  writes the rank-neutral file, then the active communicator derives ownership
  and scatters the restored hierarchy. Communicator size and owner maps are not
  serialized.
- Selected frontends permit checkpoint/stop/restart only after lexical
  input/output/checkpoint/restart collision checks. Molecular transport remains
  rejected and is not added to the schema.

## Consequences

A selected full-H2/O2 first-step checkpoint written in serial or on two ranks
is byte-identical. It resumes to the exact uninterrupted coarse and fine CSV
bytes in serial and on one, two, four, and eight ranks, including ranks that own
no fine planes. Fixed schema-2 checkpoint and output hashes remain unchanged,
and fixed/selected readers reject the other schema transactionally.

The bundle digest is provenance, not payload authentication. Formatted
root-only `status="replace"` I/O is neither crash-atomic nor scalable, and
schema migration, AMR molecular transport, dynamic 3D topology, physical
coarse boundaries, and 3D EB AMR remain separate work.
