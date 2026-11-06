# Collective native adjacent MOL donor admission, 2026-10-13

## Failure and repaired behavior

The root material-replay fixture previously timed out after 1200.24 seconds in
its first distributed weight-zero MPI case. Its serial counterpart had already
passed both unchanged gas/species and actual UEINT inventory gates. The separate
first-fine receipt replay completed all twelve cases in the preceding source;
that earlier result is not counted as a rerun of this candidate.

A bounded debugger reproduction isolates mismatched collective ordering:
`fund_exposed_coarse_endpoint` evaluates same-patch donors only on the coarse
owner, yet its native packet callback calls the collective reactor-admission
predicate. The owner waits in that predicate's `MPI_Allreduce`; other ranks have
advanced to the subsequent `MPI_Bcast`. The debugger interruption is diagnostic
evidence, not a passing numerical test.

Native adjacent admission now broadcasts each read-only target/donor proposal
and temperature/debt anchor, then evaluates the packet query on every rank in
the same order. Face overflow admission is broadcast before selecting another
donor. Accepted state, temperature, gas/scalar faces and changed flags remain
owner-local. The existing owner-only path is retained when no native callback
is supplied. Serial numerical kernels, mechanisms, reaction/EOS tolerances,
`3e-14` inventory gates and golden data are unchanged.

## Verification

A new fast `unit_mpi_eb_native_adjacent_query_3d_3` exercises a real material debt
with distributed and owner-zero layouts. Every owner/idle rank must perform the
same nonzero number of packet queries. Exact serial/MPI gas/T/UEINT and both
face families match; the existing independent signed-inventory/debt check is
retained. Rejecting the private query only on idle rank 1 collectively refuses
and retains the caller's original gas/T/UEINT and both face families bitwise.
The repaired actual full-mechanism root weight-zero reproducer also passes.
Its independent gas/species and UEINT residuals are `1.583e-17` and `1.132e-17`.

The complete root/fine material replay matrix checks twelve cases per window:
one/two MOL corrections, `[0,2]`/`[0,0]` ownership and weights 1/0.5/0. External
and auxiliary callbacks refuse; actual archived stage packets must stay frozen.
Unit weight matches original native gas/T/UEINT/reaction-source endpoints;
all weights check bitwise serial/MPI states, both face families, budgets,
funding stock/impulses and pending debts. Independent inventories include
outstanding actual receipts once at the window's entry and exit. These modes
were added in the preceding change; their full source-frozen material scope is
qualified only by the completed results below.

The final candidate passes seven CTest cases:

- `unit_mpi_eb_native_mol_material_replay_3d_3`: 769.13 seconds.
- `unit_mpi_eb_native_mol_receipt_replay_3d_3`: 399.90 seconds.
- `unit_mpi_eb_tree_flux_3d_3`: 22.91 seconds.
- `unit_mpi_eb_native_mol_replay_3d_3`: 107.81 seconds.
- `unit_mpi_eb_native_mol_weighted_replay_3d_3`: 305.86 seconds.
- `unit_mpi_eb_native_mol_subparent_replay_3d_3`: 255.05 seconds.
- `unit_mpi_eb_native_adjacent_query_3d_3`: 0.68 seconds.

All 24 full-mechanism material replay cases pass. Maximum independent
gas/species and actual UEINT inventory residuals are `5.657e-17` and
`6.379e-17`, below the unchanged `3e-14` gate. The existing weighted
and sub-parent controls retain their 12 and 24 parameter cases, including
their atomic provider/clock/policy rollback checks. Mechanism provenance,
Python provenance (14/14), syntax (186 files), project contract and
`git diff --check` pass. The preceding timeout and debugger interruption
are not counted as passes. These are focused library checks.


## Source identity

Dirty `main` remains at `0a0e00ebe9138145fbb0e0e745665d0d768df009`; existing user
changes and unrelated calculations are preserved. The GNU Fortran 13.3.0 /
OpenMPI 4.1.6 Debug configuration uses matching `/usr/bin` wrappers, enables
MPI/spray and fixed/selected `h2o2_full`, and leaves SUNDIALS OFF. The pinned PeleC
reference remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Commands, full source/linked identities, before snapshots, incremental diff,
build/debugger/green logs, JUnit and qualification are retained in
`build/native-mol-material-zero-20261013`. The original timeout, twelve-case
first-fine pass and aborted finalizer are retained separately in
`build/native-mol-material-replay-20261012`. No timed-out or interrupted run is
counted as passed. The initial query result reader found no relative JUnit
artifact; the subsecond unit was rerun with an absolute JUnit output path.
The original log and postprocessing error are retained; only the final complete
JUnit result is counted. A `zero-probe` selector narrows the real fixture for
debugger reproduction; regular CTests retain all their parameter cases.

## Remaining boundary

This qualifies explicit source-frozen root/fine replay with real material
receipts and collective adjacent-query refusal. It does not connect automatic
recovery after a partially observed native parent/future-fine failure. That
connection must retain actual observed source prefixes and must not invent a
missing future packet or retry an admitted caller-owned reactor refusal.
Prepared/routed consumers, multilevel material debts, transport/CVODE matrices,
application/restart/install audits, pinned numerical trajectories and physical
validation remain separate.
