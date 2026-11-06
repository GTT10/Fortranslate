# Strict trace-debt admission for gas/held EOS donors, 2026-10-12

## Change and evidence boundary

Both complete gas EOS and independent held EOS donor kernels now require
strictly positive density above the existing floor and nonnegative species
densities before treating a pending endpoint as already funded. A caloric EOS
query may admit a tiny negative composition within its normalization tolerance;
that does not discharge the signed conserved receipt. The unchanged search
instead selects an actual whole-material donor packet. No negative density is
clipped and no chemistry, EOS or inventory tolerance is relaxed.

The unit fixture starts with a zero H2 component and a `-1e-21` signed H2
receipt. It independently demonstrates caloric EOS admission of that endpoint,
then requires real donor transfer and strict nonnegative output in gas-only,
independent UEINT and joint gas/UEINT/reactor-energy paths. Unequal target/donor
volumes check every component against both opposite packet impulses and the
combined inventory. A donor with no H2 refuses and preserves all output arrays.
The new unit fails against the previous library, and passes after the guard fix.

The preserved original future-fine fixture is freshly linked without changing
its source, input, native builtin chemistry or `3e-14` inventory gates. Before
this change it refuses in the donor interface repair with a coarse N2 endpoint
of `-1.6653345369377349e-21`; the gas EOS kernel exits at its old caloric-only
shortcut. After this change the same interface repair reports acceptance, but
the containing root clock still refuses later. This is one corrected cause,
not qualification of the complete future-fine producer or automatic recovery.
The failing fixture and debugging records remain separate from passing tests.

## Frozen candidate

The checkout remains dirty `main`, HEAD
`0a0e00ebe9138145fbb0e0e745665d0d768df009`; user changes are preserved.
The compatible GNU Fortran 13.3.0 / OpenMPI 4.1.6 Debug build uses MPI/spray,
fixed/selected `h2o2_full`, and no SUNDIALS. The pinned PeleC reference is
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`. These focused numerical tests do
not establish installed-binary, full matrix, reference trajectory or physical
qualification.

Artifacts live in `build/native-mol-future-recovery-20261012`: before snapshots,
`incremental-v2.diff`, red/green unit logs/results, source/library identities,
unchanged future fixture source, its before/after results, and debugger logs.
An initial provenance command omitted required `--bundle` and failed with
usage status 2; the corrected byte-pinned bundle preflight passes. That command
error is retained and is not counted as a passing check.

The final trace-guard candidate passes eight CTest cases:

- `unit_eb_reactive_reflux_3d`: 0.66 seconds.
- `unit_mpi_eb_tree_flux_3d_3`: 22.27 seconds.
- `unit_mpi_eb_tree_mol_3d_3`: 1.95 seconds.
- `unit_mpi_eb_tree_mol_reflux_3d_3`: 2.04 seconds.
- `unit_mpi_eb_native_mol_replay_3d_3`: 112.56 seconds.
- `unit_mpi_eb_native_mol_weighted_replay_3d_3`: 310.91 seconds.
- `unit_mpi_eb_native_mol_subparent_replay_3d_3`: 266.68 seconds.
- `unit_mpi_eb_held_reset_ledger_3d_3`: 13.32 seconds.

The sub-parent and weighted-root controls retain their 24 and 12 parameter
cases. Mechanism provenance, Python provenance (14/14), Python syntax
(186 files), project contract and `git diff --check` pass. The unchanged future
fixture still fails in this candidate; it is not included in the passing count.
The separate staged immediate-funding experiment is not part of this source
identity and has not been promoted to the checkout.
