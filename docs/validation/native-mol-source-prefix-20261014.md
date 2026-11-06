# Observed native MOL source-prefix capture, 2026-10-14

## Completed behavior

Selected ordinary serial/MPI native MOL common steps can optionally return an
`eb_native_mol_source_prefix_3d` after a native numerical failure. This separate
type contains the actual observed external and gas/actual-UEINT auxiliary
packets, native source-stage count, interval, force/transport/chemistry flags,
policy/frame, selected entry gas/T/H/IR and available current funding stock.
It is not an accepted interval and cannot be supplied as a complete replay
history. Unobserved packet fields remain unallocated.

The optional `native_source_prefix`/`native_source_prefix_ready` pair is off by
default. MPI checks presence collectively before evaluating sources; a rank
mismatch refuses on owners and idle ranks. The initial scope is reacting,
selected ordinary clocks with built-in spatial physics, independent H and a
numerical-energy ledger. Prepared/routed wrappers, caller spatial providers,
whole clock-group prefix assembly and automatic consumption remain separate.

A native refusal publishes only the observed source cache. A declined or
malformed declared packet, an independently native-admitted caller reactor
refusal, a preflight refusal, or a successful interval leaves the previously
supplied prefix unchanged and returns readiness false. No failed physical
state, face, reaction source or accepted ledger is committed with the prefix.

## Fresh evidence

Six CTests pass:

- `unit_mpi_eb_native_mol_source_prefix_3d_3`: 11.22 seconds.
- `unit_mpi_eb_native_mol_future_source_3d_3`: 222.82 seconds.
- `unit_mpi_eb_tree_held_natural_replay_3d_3`: 77.05 seconds.
- `unit_mpi_eb_tree_held_scalar_replay_3d_3`: 79.57 seconds.
- `unit_mpi_eb_tree_held_native_funding_3d_3`: 657.36 seconds.
- `unit_mpi_eb_tree_held_native_material_3d_3`: 725.55 seconds.

The new two-layout prefix test uses the pinned full H2O2 mechanism with inert
nitrogen at 1400 K, actual cut geometry and the exposed coarse cell `(3,2,2)`.
Its successful old auxiliary packet withdraws twice the local gas/independent-H
stock. The ordinary selected clock cannot fund this density deficit, so its
native predictor/spatial admission refuses before the first new-time source.
The observed old external and auxiliary packets are retained exactly. The
three-slot/two-correction profile has one observed source stage; both future
packet slots remain absent and the caller source count is exactly one.
This is an intentionally refusing cache fixture, not an accepted calculation.

Serial and both MPI owner layouts preserve caller gas/T/H/IR, gas/scalar faces and
all numerical/source ledgers. The observed packets agree exactly. A declined new-time source, including the MPI idle rank, publishes no prefix
after the old packet was already observed. Its two original source callbacks
are retained without repeating the declined call or reaching a future stage.
The reactor refusal occurs only after the native full-mechanism solve accepts
its input; the caller is invoked once, two original source stages are observed,
and no recovery prefix is published. Rank-disagreeing optional presence
refuses before any callback. A successful two-correction native step also
matches with and without capture: gas/T/H/IR, gas/gas/scalar faces and accepted ledgers
agree exactly, with the same two original auxiliary callbacks and no prefix.

The existing current-stock source and four automatic material/parent controls
pass against the newly linked libraries. Mechanism provenance, fourteen Python
provenance tests, Python syntax, project contract and `git diff --check` pass.
The initial compiler failure is preserved in `failed-build-v1` with the actual
preprocessed inputs; it is not a passed numerical check. The first numerical launch stopped on
an unassociated test procedure pointer before entering the common driver;
its source, linked artifacts and logs remain in `failed-v1`. Default reactor
arguments are now omitted explicitly; runtime checking is unchanged. A second
candidate refused before source evaluation because the selected-step test
omitted its required gas-face output. A third exposed mismatched test assertion
collectives across owner and idle ranks; neither is counted as passed. They
remain in `failed-v2` and `failed-v3`. The new per-owner and callback checks now
use local assertions that abort the MPI job on failure. The fourth
candidate reached layout-one source/rollback checks but read an absent
inactive-level face during its success comparison. Its complete evidence
remains in `failed-v4`; it is not counted as passed. Face parity now verifies
allocation, shape and exact values; layout pass markers follow an all-rank
barrier.

## Identity and boundary

Dirty `main` is `0a0e00ebe9138145fbb0e0e745665d0d768df009`. GNU Fortran 13.3.0 /
OpenMPI 4.1.6 Debug, fixed/selected `h2o2_full`, MPI/spray and SUNDIALS OFF are
unchanged. Source snapshots/hashes, incremental four-file diff, exact build and
CTest commands, linked library/binary hashes, logs/JUnit and qualification are
in `build/native-mol-source-prefix-20261014`. Existing user work and unrelated
simulations are preserved; reaction/EOS/species and inventory gates are retained.

This qualifies observed-source capture only. Rebuilding a partially observed
native clock from these packets, assembling parent/descendant prefix windows,
requesting an unobserved source exactly once at its first native visit, and
preserving pending/open funding context through that recovery remain open.
Application/restart/install audits, broader numerical configurations, pinned
solver trajectories and physical validation are separate evidence classes.
