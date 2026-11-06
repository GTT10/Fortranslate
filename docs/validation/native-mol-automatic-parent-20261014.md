# Native automatic held MOL parent recovery, 2026-10-14

## Behavior and boundary

When a completely observed ordinary held MOL parent window fails its native
closure, serial/MPI recovery now re-enters the recursive native driver from its
pre-payment clock anchor. Each private trial freezes the actual archived
external/auxiliary stages, rebuilds coupled gas/UEINT spatial faces and chemistry,
and carries the entire pending receipt/open ancestor context. The existing
zero/unit trials and twelve bisection refinements select one common spatial
weight without changing reaction/EOS or inventory acceptance gates.

An accepted trial replaces the failed window's state, temperature, actual held
energy, reaction source, heat, both face families, physical source inventories
and captured clock history. Earlier accepted sibling inventories remain counted
once. Actual old/new predictors are retained separately from reaction sources,
so subsequent sibling clocks receive the replacement parent trajectory. The
outer closure keeps the replacement exit event and does not close or reset the
already closed native transaction again. Caller-owned provider refusals remain
terminal and roll back the entire outer transaction.

This connects recovery of fully observed windows. It does not synthesize a
missing source stage after a partially observed native failure. Automatic
selection in a window carrying nonzero material pending debt, prepared/routed
consumers, broader multilevel material debts, installed applications, pinned
solver trajectories and physical validation remain separate gates.

## Completed checks

The full-mechanism natural fixture preserves its native interval `1e-8`, one
correction, inert initial composite state and covered coarse peroxide source.
The replacement must mark all seven physical clocks as replayed and select a
reduced common weight in distributed and owner-zero layouts. Source callbacks
remain exactly 14, physical counts remain 1/2/4, and there are no composite
peroxide products or physical auxiliary inventory. Independent gas/species and
actual UEINT inventories use the original endpoint, accepted face/volume terms
and measured projection/reset events at the unchanged `3e-14` threshold.
Serial/MPI gas/T/UEINT/IR/heat and both face families agree bitwise. All native
funding clocks close and pending receipt buckets are empty. Source refusal and
native-admitted reactor refusal retain caller fields, histories, faces, ledgers
and the original funding context without another trial.

The candidate passes these ten CTests:

- `unit_mpi_eb_tree_held_native_funding_3d_3`: 620.31 seconds.
- `unit_mpi_eb_tree_held_natural_replay_3d_3`: 75.01 seconds.
- `unit_mpi_eb_tree_held_scalar_replay_3d_3`: 83.27 seconds.
- `unit_mpi_eb_native_mol_material_replay_3d_3`: 822.92 seconds.
- `unit_mpi_eb_native_mol_receipt_replay_3d_3`: 417.02 seconds.
- `unit_mpi_eb_native_adjacent_query_3d_3`: 0.78 seconds.
- `unit_mpi_eb_tree_flux_3d_3`: 23.83 seconds.
- `unit_mpi_eb_native_mol_replay_3d_3`: 125.83 seconds.
- `unit_mpi_eb_native_mol_weighted_replay_3d_3`: 327.36 seconds.
- `unit_mpi_eb_native_mol_subparent_replay_3d_3`: 270.35 seconds.

The natural automatic fixture's maximum independent gas/species and actual
UEINT residuals are `7.642e-15` and `4.448e-15`. Controls retain 24 explicit
material receipt cases, 12 weighted and 24 sub-parent cases, as well as the
ordinary automatic scalar ledger, collective adjacent-query and flux checks.
Mechanism provenance, fourteen Python provenance tests, Python syntax, project
contract and `git diff --check` pass. The preceding baseline/candidate timeouts and bounded debugger interruptions
are diagnostic evidence and are not counted as passing candidate tests.

## Identity and artifacts

Dirty `main` remains at `0a0e00ebe9138145fbb0e0e745665d0d768df009`; existing user
changes and unrelated jobs are preserved. The unchanged GNU Fortran 13.3.0 /
OpenMPI 4.1.6 Debug configuration enables MPI/spray and fixed/selected
`h2o2_full`, with SUNDIALS OFF. The pinned PeleC reference remains
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Exact commands, before snapshots and baseline library/binary hashes, incremental
diffs, full candidate source/linked identities, build logs, JUnit, independent
inventory values and qualification state are in
`build/native-mol-automatic-parent-fix-20261014`. `run.py` verifies the full source
and linked identities before publication. The preceding implementation and both 1800-second timeouts are retained in
`build/native-mol-automatic-parent-20261013`. The bounded debugger shows all
three ranks inspecting an unallocated inactive scalar face at the first coarse
clock, before automatic parent replay. The local mean-flux limiter had created
an admitted scalar family for selected clocks only. Native funding traverses
inactive clocks as possible donors; its full output therefore needs zero faces
when no donor exists. Serial/MPI now complete those absent admitted faces
before copying and checking the family. Existing material donor faces are
preserved, and nonowner fields remain unallocated.

After the absent-face repair, all seven native clocks complete but a zero-weight
root trial refuses. A second bounded debugger shows a purely absent H2O2
density of `-6.6174449004242215e-32` introduced when the native spatial RHS
cancels an opposing face delta. Built-in spatial density/species have no
volume source: their accepted rates are now assembled directly from the
blended open-face divergence. Momentum, gas-energy and actual UEINT keep
their native non-face terms, and declared diffusion providers keep their
original signed volume/face representation. Unit weight is unchanged.
No species clipping or acceptance tolerance is introduced. The preceding
14.44-second refusal and its source/linked identity remain in `failed-v1`. `build-v2.log` and the no-work
build audit identify the repaired candidate.
