# Native MOL partial root recovery, 2026-10-16

## Qualified behavior

A newly added serial/MPI fixture now accepts an ordinary native MOL clock
whose original unit-weight attempt fails before observing all source stages.
The original attempt captures two of three native stage packets. One automatic
retry at spatial weight 1/2 succeeds, retains those packets exactly and never
calls the old source twice. Complete source-frozen replay of the accepted
history, requested explicitly at weight 1/2, reproduces the automatic result
exactly and calls no source provider.

This is a single-root, no-redistribution native-library fixture. It uses the
pinned full H2O2 mechanism with inert N2 at 1400 K, a 4x4x4 cut geometry,
outflow boundaries, a 1000 m/s x velocity, a geometric x-density profile and a
1e-8 s interval with two MOL corrections. The accepted stage-zero auxiliary
source removes gas and actual H only from the last x slab at a coefficient
of 1.002 per interval; the new-stage source is zero. Built-in native spatial
physics competes with that source. There are no fine donors or off-level
material receipts, and the source is a controlled numerical diagnostic.
No production operator, EOS/reaction/conservation tolerance or golden data
changed in this qualification.

## Fresh evidence

`unit_mpi_eb_native_mol_partial_recovery_3d_3` passed on three ranks in
13.81 seconds, with zero failures and zero skips. It verifies:

- Native unit-weight refusal before the complete three-stage history, with
  absent future packets and unchanged caller gas/T/H/I_R, faces and ledgers.
- Automatic acceptance with exactly one retry, weight 1/2 and one accepted
  physical clock. The source callback is called twice in total, once for the
  old packet, and the accepted record is marked replayed.
- Exact automatic/complete-frozen-replay and serial/MPI gas/T/H/I_R,
  gas/scalar face and physical ledger agreement for owner ranks 0 and 2,
  with the other two ranks idle.
- Independent gas/species/total-energy and actual-H inventories, including
  aperture-weighted face impulses, the actual auxiliary source, H spatial
  volume work and recorded H projection/reset. Normalized residuals are
  `6.071532165918825e-18` (gas) and `5.068645156369733e-18` (H), below the
  unchanged `3e-14` criterion. The gas mass and H auxiliary impulses are
  `-4.008000000000002e-9` and `-3.333694269532385e-3` respectively.
- A valid settled funding context, no open root clock, no pending receipt
  and zero off-level funding impulse.
- Source refusal, including idle-rank refusal, and a caller reactor refusal
  after independently checking native admission. Failed outer transactions
  retain physical fields, numerical ledgers, accepted step counts and
  previously accepted nonempty faces/history/funding context.

The mechanism provenance preflight, 14 Python provenance tests, in-memory
syntax compilation of 186 Python sources, project contract and whitespace
checks also pass. The preceding eight retry/cache/current-stock/parent-window
CTests were not rerun here. Their two libraries and two test executables remain
byte-identical to the preceding successful qualification; they remain separate
dated regression evidence.

## Source identity and reproduction

Dirty `main` HEAD is `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
The qualified source manifest contains 2,883 files; the only changes from the
preceding delivery are the new Fortran test and its CMake registration.
GNU Fortran 13.3.0 / OpenMPI 4.1.6 Debug, fixed/selected h2o2_full, MPI/spray
and SUNDIALS OFF are unchanged. The new executable SHA-256 is
`7819ba5e7b9128a33c9b54919c8363af6665fcec0f5ad800b4de963a47fe2b44`.

Exact source and linked-file SHA-256 manifests, JUnit, logs and machine-readable
qualification are in
`build/native-mol-gradient-prefix-20261016/qualification`. The passing run is
`partial-root-v2.xml` / `partial-root-v2.log`. Reproduce with:

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C cmake --build \
  build/native-mol-reflux-system-clean \
  --target test_mpi_eb_native_mol_partial_recovery_3d -j2
env -i PATH=/usr/bin:/bin LC_ALL=C \
  OMPI_MCA_rmaps_base_oversubscribe=1 HWLOC_COMPONENTS=-gl,-x86 \
  ctest --test-dir build/native-mol-reflux-system-clean -V \
  -R '^unit_mpi_eb_native_mol_partial_recovery_3d_3$' --output-on-failure
```

The initial test failed to compile because its callback name collided with an
MPI symbol; that candidate is retained in `failed-build-v1`. The first compiled
test then failed an exact funding-geometry check: it compared the pre-extraction
input geometry with the tree's extracted node geometry. `failed-test-v1`
retains that source and failure. The fixture now uses the initialized node
geometry consistently for serial, MPI and funding checks. Neither failure is
counted as passed. Earlier standalone probe failures and the bounded successful
search remain exploratory evidence outside this qualification directory.

## Remaining boundary

This closes the successful partial **single common root clock** gate.
Partially observed parent/descendant windows, accepted partial recovery with
nonzero material funding, broader source models, application/installed-binary
audits, pinned-solver trajectory parity and physical validation remain open.
The full mechanism is loaded, but inert N2 does not qualify reacting mixtures.
