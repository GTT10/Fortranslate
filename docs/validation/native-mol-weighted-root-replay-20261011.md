# Native held MOL weighted root replay, 2026-10-11

## Implemented scope

`native_spatial_weight` extends the source-frozen ordinary held MOL root replay
in the serial and MPI recursive APIs. It requires a complete
`native_replay_history` and `redistribution="none"`, chemistry and paired gas /
actual UEINT. The common APIs require a source-frozen native replay interval.
One agreed finite weight in `[0,1]` applies to every physical level clock.
Prepared redistribution, routed source and Godunov replay weights stay separate.

Each fresh clock builds one first-order PCM/Lax-Friedrichs gas / UEINT face
family from its actual old state and fresh ancestor endpoints. Each native
spatial stage mixes its complete high face family with that clock's low family.
An aperture-weighted divergence of the face difference adjusts the gas and
independent scalar rates with the same weight. Native non-face terms remain in
the spatial operator; frozen external force/work, auxiliary and classified
volume source packets are retained. Weight zero copies the low faces exactly.
Weight one bypasses low construction and rate correction entirely.

This low family has inviscid faces: a nonunit weight also replaces the face part
of any native transport/provider contribution. Native wall/other non-face
contributions remain. It is an explicit conservative trial family, not a
claim of a different validated physical transport model. The caller chooses
the weight; automatic root refusal recovery is not connected by this change.

Paired native funding/admission still runs on the chosen family. A rejected
nonunit mean aborts rather than selecting a second local flux weight. The unit
case retains its original native safeguards.

The whole root forest, funding context, chemistry, reflux and reset ledgers
stay private until closure. A source/chemistry provider refusal aborts the whole
trial. Rank-disagreeing argument presence or values, NaN and out-of-range
weights reject before physics. Open sub-parent windows, partial source prefixes,
prepared/routed consumers, actual future-fine producer recovery, pinned PeleC
trajectory comparison and physical qualification remain separate work.

## Source and artifacts

The existing dirty `main` checkout is preserved at
`0a0e00ebe9138145fbb0e0e745665d0d768df009`, containing stored `origin/main`
`5a2c818899f7213dad2fbb3c6808558711aced10`. Source snapshots, incremental diffs,
source/library/binary/cache identities, build logs and execution results are
in `build/native-mol-weighted-replay-20261011`. The compatible existing build
uses GNU Fortran 13.3.0 / OpenMPI 4.1.6, Debug, fixed and selected `h2o2_full`,
MPI and spray enabled, SUNDIALS disabled. No compiler/backend configuration changed.

## Verification

The weighted fixture uses nonuniform density/velocity/temperature, pure N2
with the full H2/O2 mechanism, native forced chemistry, cut geometry, body
force and auxiliary gas energy (the independent auxiliary UEINT row is zero).
It exercises weights 1, 0.5 and 0, one/two MOL corrections, seven physical
clocks and both distributed `[0,2,0]` and owner-zero layouts with idle ranks.
It checks field/face/budget identity with serial, full fresh ancestor caches,
frozen forcing stages, changed gas and scalar face families for nonunit
weights, collective weight rejection and rollback after completed root
chemistry followed by one refusing child.

The independent composite gas/species ledger subtracts accepted face/source
and numerical impulses. The actual UEINT ledger uses its aperture-weighted
faces, independent auxiliary scalar source, local projection/closure and
wrapper reset impulses. Reacting projection already includes entry EOS and
gas-forced thermal changes; these must not be subtracted as additional gas
energy sources. Both keep the `3e-14` acceptance gate.

The final candidate passes five CTest cases:

- `unit_mpi_eb_native_mol_weighted_replay_3d_3`: 305.00 seconds.
- `unit_mpi_eb_tree_mol_3d_3`: 2.10 seconds.
- `unit_mpi_eb_tree_mol_reflux_3d_3`: 2.06 seconds.
- `unit_mpi_eb_native_mol_replay_3d_3`: 108.53 seconds.
- `unit_mpi_eb_held_reset_ledger_3d_3`: 12.74 seconds.

The weighted test completes all 12 parameter cases. Maximum independent gas
and actual UEINT residuals are `1.126e-16` and
`7.933e-17`, respectively, below the unchanged `3e-14`
gate. Unit weight matches the original native fields/faces/budgets bitwise.
A separately linked, unchanged external-stage fixture also passes both
one/two-correction cases (28.74 seconds).
Mechanism provenance, Python provenance (14/14), syntax (153 files), project
contract and `git diff --check` pass.

The first weighted attempt used an incorrect scalar diagnostic that subtracted
gas-forced energy twice; its unit-weight fields were already bitwise native.
The next attempt passed that corrected unit case and exposed the mandatory
old-UEINT halo output in low-face construction. Both failed executions are
retained and excluded from final qualification. Final source/library/test
identities are in `linked-identity-v4.json`; complete verbose logs and JUnit
results are `weighted-v4.log/xml` and `controls-v4.log/xml`. The separate cache
regression identity/result and `qualification.json` complete the record.

These are focused inviscid library checks with an inert composition. Molecular
transport, optional CVODE, reacting mixtures, restart/application consumers,
installed binaries, full configuration matrices and physical validation were
not qualified here. Historical [unit-weight evidence](native-mol-root-replay-20261011.md)
has a separate source identity.

Reproduce from the repository root using the recorded build configuration:

```sh
cmake --build build/native-mol-reflux-system-clean --target \
  test_mpi_eb_native_mol_replay_3d test_mpi_eb_tree_full_3d \
  test_mpi_eb_tree_hydro_3d test_mpi_eb_held_reset_ledger_3d \
  pelef_selected_reactive_eb_3d_runtime -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH -u DISPLAY \
  HWLOC_COMPONENTS=-gl,-x86 OMPI_MCA_hwloc_base_binding_policy=none \
  ctest --test-dir build/native-mol-reflux-system-clean -V \
  -R '^(unit_mpi_eb_native_mol_weighted_replay_3d_3|unit_mpi_eb_native_mol_replay_3d_3|unit_mpi_eb_tree_mol_3d_3|unit_mpi_eb_tree_mol_reflux_3d_3|unit_mpi_eb_held_reset_ledger_3d_3)$'
```

The archived runners intentionally hash-gate their original tested source and
linked objects. For a later rebuild, create a fresh source/binary identity and
relink the unchanged external fixture rather than bypassing those guards.
