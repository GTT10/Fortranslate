# Native held MOL sub-parent window replay, 2026-10-12

## Implemented scope

The serial and MPI recursive APIs now accept `native_replay_history` for a
complete ordinary held MOL window rooted below level zero, with an optional
`native_spatial_weight` in `[0,1]`. `dt` is the selected clock's interval;
its depth and native start time come from the first record. Root windows still
begin at time zero. Every selected/descendant physical clock must have its
complete source tables; this does not synthesize unobserved future source stages.

Recursive histories retain actual `clock_start` / `clock_interval` vectors.
Replay restores the full clock-entry gas/T/UEINT/IR/heat and numerical prefix,
all ancestor gas/T/UEINT/source endpoints, and those exact native clocks.
It does not invert interpolation fractions to guess ancestor periods. Frozen
ancestor periods, endpoints and source packets outside the selected window
remain consistent across its records. Whole-window checks validate clock
coverage, refinement intervals, interpolation fractions, owner-local fields,
clock-vector agreement, source/policy identity and the open funding context
before invoking physics.

Sub-parent windows require the existing `source_funding_context`. Its archived
pending receipts and every open gas/UEINT ancestor face family are cloned
privately before target receipt payment. The selected subtree closes through
the native driver; external ancestor clocks stay open in the returned context.
`integrated_flux` returns this window's impulses, while the returned context
already includes its contribution to the open parent family. A caller must use
one authoritative parent family rather than add the same contribution twice.
Only a complete accepting transaction publishes state, faces, histories,
ledgers and context. Metadata or native-admitted provider refusal preserves
all accepted caller outputs.

## Source and artifact identity

The preserved dirty `main` checkout is
`0a0e00ebe9138145fbb0e0e745665d0d768df009`, containing stored `origin/main`
`5a2c818899f7213dad2fbb3c6808558711aced10`. Artifacts and before snapshots are
in `build/native-mol-subparent-replay-20261012`. The compatible Debug build uses
GNU Fortran 13.3.0 / OpenMPI 4.1.6, fixed/selected `h2o2_full`, MPI and spray
enabled, SUNDIALS disabled. Compiler/backend configuration is unchanged.

The pinned PeleC reference remains
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`. This record is focused native library
evidence; it does not establish pinned reference trajectory agreement.

## Verification scope

The new three-rank fixture uses nonuniform density/velocity/temperature, pure N2
with the full H2/O2 mechanism and native forced chemistry, planar cut geometry,
body forcing and auxiliary gas energy (the independent auxiliary UEINT row is
zero). It covers one/two corrections, distributed `[0,2,0]` and owner-zero
layouts, and weights 1, 0.5 and 0 for two observed windows:

- Depth 1 at time 0: one selected clock plus its two depth-2 children.
- Depth 2 at time `1e-8 s`: one leaf clock with two ancestors still open,
  after an earlier sibling subtree has already completed.

For unit weight, the entire gas/T/UEINT/IR/heat/prefix and pending/open face
context match the next original native clock entry bitwise. Nonunit weights
change both gas and UEINT face families. Every accepted result and its fresh
ancestor payloads/clock vectors agree bitwise between serial and MPI. Fabricated
completed caller tails are replaced by the authoritative window entry.

Independent composite gas/species and actual UEINT inventories use accepted
aperture-weighted face impulses and volume/reset/projection terms at the
unchanged `3e-14` gate. For an open sub-window, these face impulses also include
exchange across its still-open coarse/fine boundary; they are not a sealed
whole-domain boundary ledger. The nominal pending receipt values are zero in
this fixture; nonzero parent face families are exercised. Nonzero material
pending debts require separate qualification.

Missing/idle-disagreeing clock vectors, altered ancestor periods, remote NaN
UEINT and partial open-face allocation reject collectively before physics.
A depth-1 provider refuses once after accepted selected-clock chemistry; a
leaf provider refuses once after its first accepted cell. Serial and both MPI
layouts retain complete fields, faces, budgets and the parent context
without retry.

The final candidate passes six CTest cases:

- `unit_mpi_eb_native_mol_subparent_replay_3d_3`: 271.77 seconds.
- `unit_mpi_eb_tree_mol_3d_3`: 2.02 seconds.
- `unit_mpi_eb_tree_mol_reflux_3d_3`: 2.05 seconds.
- `unit_mpi_eb_native_mol_replay_3d_3`: 114.56 seconds.
- `unit_mpi_eb_native_mol_weighted_replay_3d_3`: 317.21 seconds.
- `unit_mpi_eb_held_reset_ledger_3d_3`: 13.31 seconds.

All 24 sub-parent parameter cases pass. Maximum independent gas and actual
UEINT residuals are `1.239e-16` and
`9.103e-17`, below the unchanged `3e-14` gate.
The weighted-root regression completes its 12 parameter cases. The unchanged,
freshly linked external-stage fixture passes both one/two-correction cases
(27.55 seconds). Mechanism provenance,
Python provenance (14/14), syntax (153 files), project contract and
`git diff --check` pass. No interrupted or historical run is counted here.

Source/library/test hashes are in `linked-identity-v1.json`; complete verbose
logs/JUnit results are `subparent-v1.log/xml` and `controls-v1.log/xml`.
The separate cache fixture has its own source/binary identity and result.
`qualification.json`, the tested source manifest, before snapshots and
incremental diffs identify the candidate independently of older records.


## Remaining boundaries

Automatic parent/future-fine failure recovery, partial/unobserved source-prefix
production, prepared redistribution and routed consumers remain separate work.
Molecular transport, optional CVODE, reacting mixtures, nonzero pending material
debts, restart/application/installed-binary audits, full configuration matrices,
pinned solver trajectories and physical validation were not qualified here.
Earlier [root replay](native-mol-root-replay-20261011.md) and
[weighted-root replay](native-mol-weighted-root-replay-20261011.md) records retain
their own source identities.

Reproduce the focused tests from the repository root with the recorded build:

```sh
cmake --build build/native-mol-reflux-system-clean --target \
  test_mpi_eb_native_mol_replay_3d test_mpi_eb_tree_full_3d \
  test_mpi_eb_tree_hydro_3d test_mpi_eb_held_reset_ledger_3d \
  pelef_selected_reactive_eb_3d_runtime -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH -u DISPLAY \
  HWLOC_COMPONENTS=-gl,-x86 OMPI_MCA_hwloc_base_binding_policy=none \
  ctest --test-dir build/native-mol-reflux-system-clean -V \
  -R '^(unit_mpi_eb_native_mol_subparent_replay_3d_3|unit_mpi_eb_native_mol_weighted_replay_3d_3|unit_mpi_eb_native_mol_replay_3d_3|unit_mpi_eb_tree_mol_3d_3|unit_mpi_eb_tree_mol_reflux_3d_3|unit_mpi_eb_held_reset_ledger_3d_3)$'
```
