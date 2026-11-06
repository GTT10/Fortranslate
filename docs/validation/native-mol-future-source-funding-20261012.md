# Native MOL current-stock funding and actual future fine entry, 2026-10-12

## Implemented behavior

Ordinary native gas/actual-UEINT admission with the supplied private projection
query now budgets the represented current physical states and their read-only
pending source receipts. It does not spend an unobserved future coarse/fine
reflux credit before reaction. The serial and MPI paired funding kernels route
actual conserved whole-material packets, including an adjacent covered coarse
proxy's physical descendant inventory, and retain both copies of their gas/H
face impulses. Current-clock source debt and deferred donor receipts are charged
once. The existing closure reconciliation path remains selected when there is
no native projection query. MPI includes query presence in its rank-uniform
metadata contract before choosing the transaction path.

The prior [trace-debt guard](native-mol-trace-debt-admission-20261012.md) removes
an earlier caloric-only shortcut. After that fix, interface repair accepted but
the actual root delivery still had N2 density `-0.15133433751678527`. The frozen
staged experiment separates the next cause: adding only the serial current-stock
change passes serial gas/H closure but leaves MPI refusing; the matching MPI
change passes both one/two-correction paths with exact owner outputs. These
records are diagnostic evidence, separate from the final repository tests.

## Fixture and acceptance criteria

`test_mpi_eb_native_mol_future_source_3d` keeps the preserved fixture's full
H2/O2 mechanism, builtin forced chemistry, 1400-K material packet, planar cut
geometry and source schedule. A two-level tree starts with nitrogen-free exposed
coarse neighbors and density-4 fine donors. On the first corrector, a physical
root cell withdraws 0.9 density of a pure-N2 gas/actual-UEINT packet. The first
native fine clock must restore that root-produced nonzero gas/H receipt from
its actual pre-entry anchor, pay it exactly once before building the predictor,
and remove the paid active debt from read-only stock. The next clock observes
no duplicate receipt. One/two MOL corrections and `[0,2]` / `[0,0]` owner layouts
exercise one/two idle ranks. The later Picard policy retains the zero old
auxiliary source in the accepted two-correction RHS.

Complete native serial/MPI gas/T/UEINT, gas/scalar faces, physical source and
chemical budgets, numerical prefixes, archived stocks/receipts and clock counts
agree bitwise. Every archived pre-entry context is valid at its clock depth;
the completed root returns no open family or pending stock. Independent
composite gas/species and actual UEINT inventories retain the `3e-14` gate.
The test source is extracted from the unchanged preserved fixture, with unused
historical checks removed and the second owner layout and exact old-UEINT
payment assertion added. No mechanism, reaction tolerance, source magnitude,
donor density or inventory threshold is adjusted to gain acceptance.

## Source and build identity

Dirty `main` remains at `0a0e00ebe9138145fbb0e0e745665d0d768df009`; existing
user changes and unrelated jobs are preserved. The compatible GNU Fortran
13.3.0 / OpenMPI 4.1.6 Debug build uses MPI/spray and fixed/selected `h2o2_full`,
with SUNDIALS disabled. Compiler/backend settings are unchanged. The reference
is still PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Artifacts are in `build/native-mol-future-recovery-20261012`. Before snapshots,
staged serial/MPI diffs and their distinct failure/pass records, final source
and linked identities, complete logs/JUnit, and the preserved original fixture
make the result reproducible independently of the prior root/sub-parent records.
The scalar-payment assertion is the only test-only difference between the v3
and v4 source identities; the production code is identical.

The final current-stock funding candidate passes eleven CTest cases:

- `unit_mpi_eb_native_mol_future_source_3d_3`: 214.78 seconds.
- `unit_eb_reactive_reflux_3d`: 0.68 seconds.
- `unit_mpi_eb_tree_flux_3d_3`: 24.19 seconds.
- `unit_mpi_eb_tree_mol_3d_3`: 2.10 seconds.
- `unit_mpi_eb_tree_mol_reflux_3d_3`: 2.07 seconds.
- `unit_mpi_eb_tree_held_auxiliary_3d_3`: 24.79 seconds.
- `unit_mpi_eb_tree_held_natural_replay_3d_3`: 69.54 seconds.
- `unit_mpi_eb_native_mol_replay_3d_3`: 111.63 seconds.
- `unit_mpi_eb_native_mol_weighted_replay_3d_3`: 312.74 seconds.
- `unit_mpi_eb_native_mol_subparent_replay_3d_3`: 266.04 seconds.
- `unit_mpi_eb_held_reset_ledger_3d_3`: 12.66 seconds.

All four future-fine parameter cases pass, including two nonzero material
receipt entry cases. Maximum independent gas and actual UEINT residuals are
`8.782e-18` and `4.405e-18`, below the unchanged `3e-14` gate. The
sub-parent and weighted-root controls retain their 24 and 12 parameter cases.
Mechanism provenance, Python provenance (14/14), syntax (186 files), project
contract and `git diff --check` pass. Staged runs and older eight-case results
are not counted as final-source passes.


## Boundaries

This closes the observed root source-admission refusal and qualifies actual
native future fine receipt payment in this fixture. It does not connect
source-frozen weighted whole-parent recovery after a partially observed native
failure, nor qualify frozen replay with nonzero material pending receipts.
Prepared/routed consumers, molecular transport, optional CVODE, multilevel
future debts, restart/application/installed-binary audits, complete configuration
matrices, pinned reference trajectories and physical validation remain separate.
