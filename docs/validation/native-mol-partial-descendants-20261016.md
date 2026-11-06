# Native MOL partial root recovery with descendants, 2026-10-16

## Qualified behavior

The ordinary native MOL recursive driver now has a two-level regression for
successful recovery of a partially observed root clock followed by both child
clocks. The original unit-weight root attempt observes two of three native
source stages before refusing. One automatic common-clock retry at spatial
weight 1/2 succeeds. Its child then completes two ordinary half intervals,
retaining the still-open root funding/face context until synchronization.
Only the root record is replayed; the two child records are not replayed.

This extends the preceding
[single-root qualification](native-mol-partial-root-recovery-20261016.md).
The pinned full H2O2 mechanism is loaded with inert N2 at 1400 K. Root density
follows `2**(3-i)` along x, velocity is 1000 m/s, the root interval is 1e-8 s
and both levels use two MOL corrections with native spatial physics. A child
refines coarse x cells 1:2 and the full y/z extent by two. The old auxiliary
packet withdraws gas and actual H only in root x slab 4 at coefficient 1.002
per interval; all subsequent and child packets are zero. Boundary faces are
outflow. Independent H, ordinary redistribution and the existing strict
energy policy are unchanged.

No runtime numerical source, mechanism, tolerance or golden data changed.
The source changes for this gate are the new test and its CMake registration.

## Fresh checks

`unit_mpi_eb_native_mol_partial_descendants_3d_3` passed in 37.85 seconds on
three ranks, with zero failures and zero skips. It verifies:

- The original common root clock refuses with a two-stage prefix, an absent
  future packet and unchanged caller fields, faces, receipts and ledgers.
- Automatic success uses exactly one retry and six source calls over the
  whole forest, including exactly one old-root call. Accepted physical counts
  are one root and two child clocks.
- Distributed owners `[0,2]` and owner-zero `[0,0]` reproduce serial gas,
  temperature, actual H, reaction source, both physical face families,
  funding packets and physical/numerical ledgers bitwise. Idle fields retain
  their sparse allocation schema.
- Both child pre-entry anchors preserve an open root face family. All final
  funding buckets are zero and every face clock is closed. This fixture has
  exactly zero off-level material receipts.
- Independent composite gas/species/total-energy and actual-H inventories
  use exposed leaf volume, aperture-weighted accepted faces, the original
  physical source, each native local H volume rate and recorded root
  projection/reset events. Normalized residuals are
  `6.765421556309548e-17` and `8.619407271259760e-18`, below the existing
  `3e-14` gate. Reported boundary and source budgets are also independently
  reconstructed. The gas mass/H withdrawals are `-4.008000000000002e-9`
  and `-3.333694269532385e-3`.
- A source refusal during the last child's new stage, including a refusal
  exclusively on idle rank 1, rolls back the entire root transaction after
  its root recovery and first completed child. A reactor refusal at the same
  late stage follows an independent successful native solve and is invoked
  exactly once on the child owner. Both failures retain caller gas/T/H/I_R,
  physical/numerical ledgers, accepted counts, previous nonempty gas/scalar
  faces and all selected-level clock states, plus a valid settled context.

Mechanism provenance, 14 Python provenance tests, in-memory syntax compilation
of 186 Python sources, project contract and `git diff --check` pass. Earlier
single-root and eight retry/current-stock/parent-window CTests were not rerun:
their two libraries and three test executables remain byte-identical to the
preceding successful deliveries. Those are separate dated evidence classes.

## Identity and reproduction

Dirty `main` remains at `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
The qualified source manifest contains 2,885 files. GNU Fortran 13.3.0 /
OpenMPI 4.1.6 Debug, fixed/selected h2o2_full, MPI/spray and SUNDIALS OFF
are unchanged. The new executable SHA-256 is
`07464cda874f77422f15eb4e282710d3ad61756bd345910c687611783c62b6a3`.

Artifacts are in `build/native-mol-partial-descendants-20261016`:
source/linked SHA-256 manifests, before snapshots, build logs,
`descendants-v2.log`, `descendants-v2.xml`, runner commands/result and
`qualification.json`. Reproduce with:

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C cmake --build \
  build/native-mol-reflux-system-clean \
  --target test_mpi_eb_native_mol_partial_descendants_3d -j2
env -i PATH=/usr/bin:/bin LC_ALL=C \
  OMPI_MCA_rmaps_base_oversubscribe=1 HWLOC_COMPONENTS=-gl,-x86 \
  ctest --test-dir build/native-mol-reflux-system-clean -V \
  -R '^unit_mpi_eb_native_mol_partial_descendants_3d_3$' --output-on-failure
```

The first candidate failed its common-root preflight: the test populated the
pending signed receipt input with initial fine gas instead of the empty
pre-entry queue. The corrected test supplies zero pending receipts, matching
the recursive wrapper. `failed-v1` retains the failed source, identities and
JUnit; that failure is not counted as a pass. No production admission rule
was changed. The earlier bounded two-level search was a candidate-finding
probe, not qualification.

## Remaining boundary

This qualifies a recovered partial **root common clock** with subsequently
completed descendants, including late whole-transaction rollback. It does not
qualify a partial fine/sub-parent clock recovered under an already-open
ancestor, partial recovery with nonzero pending material receipts, or broader
partially observed window replacement. Prepared/routed sources, reacting
mixtures, application/restart/install audits, pinned-solver trajectories and
physical validation remain separate gates.
