# Automatic material sub-parent MOL recovery, 2026-10-14

## Qualified behavior

The preceding completed-window automatic recovery fixture carried no resolved
material receipt. This controlled three-level full-H2O2 fixture requires actual
current-stock nitrogen and independent UEINT borrowing before automatic
sub-parent recovery. The already qualified runtime libraries are unchanged;
this unit adds a source/receipt qualification mode and its automated CTest. A separate `refusal-probe` argument runs only the reactor
rollback path; its empty entry histories are diagnostic evidence, while the
full CTest also verifies retention of previously accepted nonempty histories.

The exposed coarse target is cell `(3,3,2)`. It starts at density `0.93`,
135000 Pa and equal nitrogen/argon mass fractions. Other exposed coarse cells
contain argon, while fine donors contain nitrogen at density 4 and the target's
consistent temperature. A new-time coarse auxiliary stage withdraws a pure
nitrogen gas/actual-UEINT packet at rate `1.8e8`; the actual `1e-8` clock and
native half-stage weight withdraw density `0.9`. The initial target nitrogen
stock is only `0.465`, so actual future fine material must fund the deficit.
The source uses the independently evaluated nitrogen gas/UEINT at the target
temperature. Actual current donor stock, not unobserved future reflux credit,
is the borrowing authority.

A separate peroxide source acts only in covered fine cell `(3,3,2)` of node 2,
which has real grandchild consumers. That virtual predictor source causes the
fine closure failure under the still-open material-funded root. The root
remains unreplayed; both fine/sub-parent windows and their four grandchild
clocks must be replaced. A reduced common spatial weight is selected without
changing any species, reaction/EOS or `3e-14` inventory gate.

## Completed checks

The new mode and three controls pass:

- `unit_mpi_eb_tree_held_native_material_3d_3`: 867.66 seconds.
- `unit_mpi_eb_tree_held_natural_replay_3d_3`: 78.13 seconds.
- `unit_mpi_eb_tree_held_scalar_replay_3d_3`: 83.20 seconds.
- `unit_mpi_eb_tree_held_native_funding_3d_3`: 791.90 seconds.

The replacement histories retain resolved nitrogen and UEINT donor debits and
a nonzero pending pre-payment anchor. The first fine clock's entry species
and independent UEINT equal their pre-entry values plus the actual root
receipt exactly once. All final pending buckets are zero and all funding
face clocks are closed. Maximum archived loan/pending magnitudes are
`[52402.534676832205, 52402.534676832205]`; the test also separately requires resolved nitrogen and UEINT
rows, not merely a nonzero energy field.

Serial/MPI states, temperature, UEINT, reaction source, heat and physical gas,
scalar and registration faces agree bitwise in distributed and owner-zero
layouts. Auxiliary source callbacks remain exactly 14 and accepted physical
counts remain 1/2/4. The accepted gas/UEINT auxiliary source matches the original
physical withdrawal independently; the covered peroxide source creates no
composite chemical products. Maximum independent composite gas/species and
actual UEINT residuals are `7.484e-15` and `6.681e-15`, below `3e-14`. The H inventory
uses original/accepted H, complete accepted faces, archived actual local
pressure/volume rates and measured projection/reset events. Its source is
already present in those local rates and is not counted a second time.

The source-refusal test uses the same default integrator and withdrawal as the
successful material transaction. The separate reactor-refusal fixture withdraws
nitrogen at `2.2e8`, giving a half-stage density withdrawal of `1.1` from the
initial `0.93`. This density deficit requires native current-stock admission
before any caller reactor is invoked. A supplied reactor owns its signed-species
solve and is not previewed to request composition-only borrowing. All root and
fine callbacks in this fixture use the native full-mechanism integrator; the
intended refusal is injected only after its fine input succeeds. The accepted
material inventory fixture remains at `1.8e8`.

A caller source refusal and an independently native-admitted reactor refusal
inside the final automatic sub-parent recovery preserve caller gas/T/UEINT/IR,
heat, histories, scalar faces, numerical/source ledgers and the original
funding context without another weighted trial. The reactor rejection occurs
after the source schedule is fully observed and after earlier native clocks
have completed. Mechanism provenance, fourteen Python provenance tests, Python
syntax, project contract and `git diff --check` pass.

## Source and remaining boundary

Dirty `main` remains at `0a0e00ebe9138145fbb0e0e745665d0d768df009`. GNU Fortran
13.3.0 / OpenMPI 4.1.6 Debug, fixed/selected `h2o2_full`, MPI/spray and SUNDIALS
OFF are unchanged. Both runtime library hashes must match the preceding
[automatic recovery candidate](native-mol-automatic-parent-20261014.md).
Existing user changes and unrelated simulations are preserved.

Exact source/linked identities, commands, source inputs, before snapshots,
incremental test diff, build/CTest logs, JUnit and qualification are in
`build/native-mol-automatic-material-20261014`. This is focused native library
evidence for one controlled three-level material topology. Partially observed
source/chemistry recovery, prepared/routed consumers, broader material debt
networks, application/restart/install audits, pinned numerical trajectories,
optional backend matrices and physical validation remain separate gates.

The initial 462-second candidate passed both owner-layout material inventories
but failed its source-refusal counter. A second 474.685-second candidate passed
that source refusal but failed to reach the intended reactor-refusal point.
The focused diagnostic observed two source calls and zero reactor callbacks.
A bounded debugger capture found the root flux-only limiter refusing before
reaction, with native material funding inactive. Those candidates remain in
`failed-v1` and `failed-v2` and are not counted as passing CTests. The refusal
fixture now requires density funding explicitly; the proposed test-only analytic
inert root adapter was removed without qualifying it.
