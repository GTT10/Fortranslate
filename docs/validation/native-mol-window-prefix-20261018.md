# Native MOL window-prefix continuation, 2026-10-18

The serial and MPI recursive ordinary held-MOL APIs accept the optional
`native_replay_prefix=.true.` together with `native_replay_history`.
Previously, the supplied records had to cover every descendant clock through
the end of the requested window. The new opt-in mode accepts a nonempty
leading sequence of complete clock records, replays their observed external
and auxiliary source packets, and evaluates subsequent clocks natively.

This is a prefix of complete clocks, distinct from the intra-clock
`eb_native_mol_source_prefix_3d` cache used by automatic common-step retries.
It does not synthesize unobserved packets or skip a missing clock in the
middle of the supplied sequence. A supplied spatial weight applies only to
replayed clocks; newly reached clocks start with their native spatial weight
and retain the existing bounded source-preserving retry path. Whole-window
parent repair is still disabled during explicit replay.

The default remains strict complete-window replay. An empty prefix, a true
prefix flag without history, or rank-disagreeing flag values is rejected.
Clock depth/time/interval are checked against the recursive schedule before
using each cached record. Gas, temperature, actual internal energy, reaction
history, physical faces, pending receipts and ledgers remain private until
the entire window succeeds.

## Focused validation

`unit_mpi_eb_native_mol_partial_subparent_3d_3` passed in 191.015 seconds,
with zero failures and zero skips. The archived source hashes were rechecked
on 2026-10-19 before further changes.
The existing three-level native recovery fixture is extended with a root-only
history, six newly reached descendant clocks, exact native/prefix final
fields/faces/ledgers, source-call counts (12 instead of 14), retry of the
first fine clock, late descendant source refusal, strict truncated-history
rejection and MPI flag disagreement. Two layouts include an idle rank.
The existing independent gas/H inventory gate remains unchanged. Mechanism
provenance, 14 Python provenance tests, in-memory Python syntax compilation,
project contract and diff whitespace checks passed in this delivery.
The initial build caught line truncation in the added Fortran statements;
those statements were wrapped before the successful build and test.

## Reproduction and scope

Dirty main base: `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
GNU Fortran 13.3.0 / OpenMPI 4.1.6, Debug. Existing build configuration is
unchanged: fixed/selected h2o2_full, MPI/spray ON, SUNDIALS OFF.
Incremental source diffs, before snapshots, build log and checks are under
`build/native-mol-window-prefix-20261018`.

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C cmake --build \
  build/native-mol-reflux-system-clean \
  --target test_mpi_eb_native_mol_partial_subparent_3d -j2
env -i PATH=/usr/bin:/bin LC_ALL=C \
  OMPI_MCA_rmaps_base_oversubscribe=1 HWLOC_COMPONENTS=-gl,-x86 \
  ctest --test-dir build/native-mol-reflux-system-clean -V \
  -R '^unit_mpi_eb_native_mol_partial_subparent_3d_3$' --output-on-failure
```

This closes the library connection from a complete-clock prefix into native
unobserved descendants. Nonzero-material partial recovery, prepared/routed
sources, arbitrary window replacement, application/restart/install checks,
pinned solver trajectories and physical validation remain separate gates.
