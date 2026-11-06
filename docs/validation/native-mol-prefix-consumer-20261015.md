# Explicit native MOL source-prefix consumption, 2026-10-15

## Behavior and boundary

Selected ordinary reacting serial/MPI MOL common steps accept an optional
`native_source_replay` of type `eb_native_mol_source_prefix_3d`. A replay consumes
only observed external and gas/actual-UEINT auxiliary packets. It rebuilds the
native spatial family on the current physical caller state and evaluates an
unobserved source only at its first reached native stage. The existing common
spatial weight can also select the native/PCM face family with a prefix input.
No accepted state, face family, reaction source or ledger is taken from the
failed physical trial.

The consumer checks the selected depth, time interval, correction count,
chemistry/transport flags, force/acceleration, frame, dual-energy policy,
topology and geometry, post-entry-reset gas/T/H/IR, correction and current funding
stock. MPI checks option presence and observed-stage count collectively, and
owner-local field validation includes every owned level. Future packet fields
must remain absent. Prefix capture now retains whole owned entry fields and
native source geometry/acceleration to make this binding possible.

A successful callback is recorded before native spatial evaluation. A later
native failure therefore retains that actual packet even when the spatial
operator refuses. Publishing a longer prefix merges previously observed input
packets and removes incomplete future slots. A caller source refusal or a
native-admitted reactor refusal publishes no new prefix. Transactional caller
state, faces and ledgers are preserved on all refusals.

This is explicit common-clock packet reuse. Recursive pre-payment anchors,
ancestor history, parent/descendant prefix assembly, pending/open funding queues
and automatic partial-window recovery are not connected by this increment.
Prepared/routed providers, wider configurations, installed binaries,
application/restart audits, pinned external solver trajectories and physical
validation are separate evidence classes.

## Focused evidence

Seven CTests passed, with zero failures and zero skips:

- `unit_mpi_eb_native_mol_source_prefix_3d_3`: 10.82 seconds.
- `unit_mpi_eb_native_mol_prefix_consumer_3d_3`: 45.96 seconds.
- `unit_mpi_eb_native_mol_future_source_3d_3`: 222.88 seconds.
- `unit_mpi_eb_tree_held_natural_replay_3d_3`: 74.81 seconds.
- `unit_mpi_eb_tree_held_scalar_replay_3d_3`: 76.61 seconds.
- `unit_mpi_eb_tree_held_native_funding_3d_3`: 615.55 seconds.
- `unit_mpi_eb_tree_held_native_material_3d_3`: 714.41 seconds.

The new consumer test uses the pinned full H2O2 native reactor with inert N2 at
1400 K, existing two-level cut geometry and owners `[0,2]` and `[0,0]` on three
MPI ranks. It checks serial and owner-local common steps. Successful reuse
controls slice the first one or two packets from an actual successful native
producer's complete source record; they do not claim a naturally recovered
partial failed trajectory. The supplied observed prefix metadata comes from the
native refusal capture fixture on the same clock entry. Nonzero source packets
are recorded by the native producer rather than fabricated as absent future
stages.

At unit spatial weight, one observed stage causes exactly one new source call,
and two observed stages cause none. Gas/T/H/IR, gas/scalar faces and all ledgers
match the fresh native reference exactly. At spatial weight 0.5, prefix reuse
matches the existing complete frozen-source weighted reference exactly.
Observed old callbacks explicitly abort if invoked again.

A new-source refusal on an idle rank, a missing observed packet, a populated
future packet and rank-mismatched optional input refuse without publication or
caller mutation. The original twice-stock sink remains a native refusing
prefix and makes no new source call. A native-admitted caller reactor refusal
also publishes no new prefix. A newly accepted sink packet followed by native
chemistry refusal extends the cache to two observed stages; consuming that
cache again makes no caller source calls and preserves the same two-stage
observation with no invented third packet.

The existing capture, current-stock funding and four automatic complete-window
controls are independent regression checks. Their success does not establish
whole-PeleF or automatic partially observed wrapper recovery.

## Provenance

Dirty `main` starts at `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
GNU Fortran 13.3.0 / OpenMPI 4.1.6 Debug, fixed/selected `h2o2_full`, MPI/spray
and SUNDIALS OFF are unchanged. Exactly four source/test/CMake files changed
from the previously qualified source-prefix delivery; no mechanism, EOS or
conservation threshold was changed. Source SHA-256, the incremental diff,
linked identities for two libraries/two executables, command/result records,
logs/JUnit and qualification live in `build/native-mol-prefix-consumer-20261015`.
The superseded pre-qualification candidate is retained separately and is not
counted as passed. Fourteen Python provenance tests, Python syntax,
mechanism provenance, project contract and diff whitespace are required.
