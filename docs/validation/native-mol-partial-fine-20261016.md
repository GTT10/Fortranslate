# Native MOL partial fine-clock recovery, 2026-10-16

## Qualified behavior

A two-level ordinary native MOL regression now recovers a partially observed
fine clock while its root ancestor's face/funding family remains open. The
root and second fine records are not replayed; only the first fine clock
uses one automatic retry at spatial weight 1/2. The complete forest retains
one root and two fine physical clocks and exactly six source callbacks,
including one old packet per fine clock.

The fixture has exposed coarse cut cells and a fine patch covering coarse
x cells 3:4 and the full y/z extent at refinement ratio two. Fine x densities
are 4, 2, 1 and 1/2; the covered coarse means are 3 and 3/4, while exposed
coarse fluid has density 8. Gas is inert N2 at 1400 K and 1000 m/s with the
pinned full H2O2 mechanism, native spatial physics, two MOL corrections and
outflow faces. The root interval is 1e-8 s. Each fine old auxiliary packet
withdraws gas and actual H only from its last x slab at coefficient 1.002
per fine interval; other packets are zero. H is initialized independently
from thermodynamics. No production operator, tolerance, mechanism or golden
data changed.

This extends the preceding
[partial root with descendants](native-mol-partial-descendants-20261016.md).
It exercises a leaf fine clock under an open parent, rather than a sub-parent
with descendants or a material-funded recovery.

## Fresh evidence

`unit_mpi_eb_native_mol_partial_fine_3d_3` passed in 58.83 seconds on three
ranks, with zero failures and zero skips. Its checks include:

- Exact serial/MPI gas, temperature, actual H, reaction source, both physical
  face families, funding packets and physical/numerical ledgers for owners
  `[0,2]` and `[0,0]`. Accepted replay flags are `[false,true,false]`, the
  retry count is one and physical counts are `[1,2]`.
- Independent reconstruction of the original first fine clock from its
  archived entry, gas/H/source ancestor endpoints, time fractions and exact
  signed funding stock. Its unit-weight native attempt refuses after two of
  three source stages, leaves the future packet absent and preserves caller
  fields, faces, receipts and numerical ledgers.
- Explicit weight-1/2 consumption of that independently captured prefix
  matches the automatic first fine clock's gas/T/H/I_R and gas/H faces
  bitwise, with no repeated source callback. Both source packets observed
  before refusal equal those in the automatically accepted history.
- Both fine pre-entry anchors retain the open root face clock. Final pending
  receipts are zero and every face clock closes. All accepted off-level
  material receipts are exactly zero in this fixture.
- Independent composite inventories over exposed leaf volume. Gas/species
  and total energy use actual accepted aperture-weighted physical faces,
  source/chemistry/wall terms and energy correction. Actual H additionally
  uses each native local volume rate and measured projection/reset events.
  Normalized residuals are `1.327063459122257e-16` (gas) and
  `5.735429492448319e-17` (H), below the unchanged `3e-14` gate. The source
  budget is independently reconstructed from both physical fine entry states;
  integrated gas mass/H withdrawals are `-3.023063666739372e-9` and
  `-2.514458108267134e-3`.
- A last-fine new-stage source refusal, including an idle-rank-only refusal,
  and a reactor refusal after an independent successful native solve roll
  back the entire transaction after the recovered first fine clock. They
  preserve all caller gas/T/H/I_R, physical/numerical ledgers, accepted counts,
  prior nonempty gas/H faces and all selected-level accepted clock states,
  plus a valid settled funding context. The late reactor is declined exactly
  once on the fine owner.

Mechanism provenance, 14 Python provenance tests, in-memory syntax compilation
of 186 Python sources, project contract and `git diff --check` pass. Previous
retry/root/descendant CTests were not rerun here. Their two libraries and four
executables remain byte-identical to the preceding deliveries; those results
remain separate dated evidence.

## Identity and reproduction

Dirty `main` HEAD remains `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
The qualified source manifest contains 2,887 files, changing only the new
Fortran test and CMake registration from the preceding delivery. GNU Fortran
13.3.0 / OpenMPI 4.1.6 Debug, fixed/selected h2o2_full, MPI/spray and SUNDIALS
OFF are unchanged. The new executable SHA-256 is
`ec1231348f5a0923eda506cbfbec740fa395c6164171ca6feaeeb76f349ea7b1`.

Exact source/linked hashes, before snapshots, commands/results, build log,
`fine-v1.log`, `fine-v1.xml` and `qualification.json` are in
`build/native-mol-partial-fine-20261016`. Reproduce with:

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C cmake --build \
  build/native-mol-reflux-system-clean \
  --target test_mpi_eb_native_mol_partial_fine_3d -j2
env -i PATH=/usr/bin:/bin LC_ALL=C \
  OMPI_MCA_rmaps_base_oversubscribe=1 HWLOC_COMPONENTS=-gl,-x86 \
  ctest --test-dir build/native-mol-reflux-system-clean -V \
  -R '^unit_mpi_eb_native_mol_partial_fine_3d_3$' --output-on-failure
```

The exploratory full-coverage fine probe initially stopped before any trial
because an inherited root-leaf assertion did not apply to that topology.
Its failed source/log/identity are preserved under the preceding delivery's
`fine-probe/failed-v1`. Corrected full-coverage and mixed-leaf probes each
found a candidate on the first trial; they remain exploratory. The promoted
mixed-leaf CTest passed on its first run.

## Remaining boundary

Partial sub-parent clocks with descendants, accepted partial recovery with
nonzero pending gas/H material receipts and broader partial-window replacement
remain open. Prepared/routed sources, reacting mixtures, application/restart/
install audits, pinned solver trajectories and physical validation remain
separate gates.
