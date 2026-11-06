# Opt-in EB reaction admission — 2026-11-06

## Contract

The serial and MPI EB library entrypoints accept an optional
`evaluate_reaction_admission` capability alongside the existing side-effecting
`evaluate_reaction` provider, plus optional `reaction_admission_status` output.
The legacy provider signature and the no-capability path remain available.
This is a library option, not an application input or a changed global solver
policy. Chemistry and an actual evaluator are required. Nonheld reaction and
held-energy common/replay/subcycle paths forward the capability.

Admission receives immutable species densities, temperature seed, initial held
energy, gas-derived external energy, species forcing, complete duration and
unchanged tolerances. The ordinary `react_sdc_state_3d` packing is used for
private EB donor queries. Gas and held energies are not equated. The callback
must be side-effect-free and repeatable for a fixed packet and fixed
species/mechanism/provider configuration during one advance.

| Status | Meaning |
| --- | --- |
| `reaction_admission_not_attempted` | No admission query was performed |
| `reaction_admission_accepted` | The supplied capability certifies the whole requested interval |
| `reaction_admission_rejected` | The supplied capability has a physical admissibility certificate |
| `reaction_admission_unknown` | No certificate; includes numerical failure or exhausted work budget |
| `reaction_admission_provider_refused` | Explicit refusal; cannot be overridden by donor recovery |

`check_forced_reactor_admission` is a copied-state builtin helper with a required
caller-selected `maximum_steps`. It returns only accepted or unknown. A false
legacy solver result does not reveal a physical certificate. Its work budget
does not change the ordinary reactor's default 100000-step limit.

Private conservative funding can try another proposal after unknown; it cannot
interpret unknown as physical impossibility or silently discard source debt.
A complete funded family is rechecked before publication. Immediately before
each actual evaluator invocation the exact assembled packet is admitted again,
without memoization. An evaluator refusal after that certification is terminal;
no substitute builtin solve or retry can override it.

Each advance has at most 16 cached private-query results. Keys use the bitwise
representations of all 2N+6 packet scalars. The table, mechanism and capability
configuration are fixed within the advance. Both accepted and unknown are
cached unchanged; allocation failure simply disables caching. This does not
change the conservative search, tolerances or reactor arithmetic. MPI agrees
presence and terminal status across owners and idle ranks.

## Conservative donor admission

The opt-in path also checks native reacting endpoints when selecting exposed
targets. A statically admissible gas/held-energy endpoint can still need funding
for the actual forced reaction trajectory. Previously that target was skipped
on a later SDC iteration. Coarse and fine previews retain the exact source-debt
mapping used by publication; the source is included once.

The shared donor helper has a second opt-in check before its statically valid
no-work return. An accepted native target keeps that return. An uncertified
target enters the existing conservative donor search. Both switches default to
false, and are enabled only with the supplied capability in the corrector.
No source, timestep, mechanism, tolerance or reactor default was changed.

## Completed focused checks

The final frozen integration passed **17/17 CTests, with no failures or skips**:
four opt-in kinetic cases (Godunov/MOL, one/two SDC iterations), the canonical
inventory regression, and the twelve contract/conservation/refusal checks.
The original exact-state, conservation and Godunov replay assertions passed.
Source and executable hashes were verified unchanged across the run.

The integration tree is based on 3073 exact files from the phase2-applied dirty
main snapshot. A fresh GNU 13.3/OpenMPI Debug CMake build produced the relevant
core and MPI archives. The initial new serial test exceeded the standard
132-column source limit; wrapping its two long lines resolved the build error.
No compiler-wide relaxation was introduced.

The twelve focused regression checks cover: immutable helper,
signed replenishment, serial provider/rollback and native donor contracts,
legacy inert material and manufactured caloric provider, opt-in inert material,
and MPI unknown/explicit refusal/actual-provider refusal/missing provider/
mismatched capability checks. MPI tests include an idle rank. Allocated flux,
state and ledger sentinels survive refusal. Independent mutation checks reject
false physical-rejection certificates, the old premature donor return and an
incorrect post-debt preview packet.

The observed original donor remains preserved as an independent exact-packet
regression. With 512 steps it returns unknown and leaves its inputs untouched.
The earlier hanging legacy kinetic invocation is not relabeled as a pass.

## Kinetic qualification status

The test-only admission policy first uses 512 steps. An unknown result with
positive external energy and positive net mass forcing is retried with 8192
steps. Other uncertified packets remain unknown. This policy applies only to
the opt-in fixture; it does not change the ordinary solver or legacy caller.

The first four-case root CTest run completed with three passes and one failure.
Godunov iterations 1 and 2 passed, including recorded-interval replay; MOL
iteration 2 passed. MOL iteration 1 accepted both serial and MPI advances and
their individual conservation ledgers, but failed the combined serial/MPI
state/energy/ledger/pressure-work comparison.

A diagnostic found bit-identical actual reactor inputs and exact gas,
temperature and held-energy arrays. Only chemical component 13 exceeded the
unchanged 32-epsilon comparison: the difference was 7.222237291452134e-35
against a 6.618668364594347e-35 bound. Serial integration accumulated cells
continuously across nodes while MPI summed node totals. An independent
arithmetic oracle reproduced all ten species ledger components under both
groupings, including the discrepant component bit-for-bit.

The serial inventory now uses the existing MPI convention: sum cells within
each node, then sum node totals in tree order. A cheap signed cancellation
regression fails the old grouping and requires exact serial/MPI agreement,
covered-cell exclusion, swapped ownership and an idle rank. The final integrated
run passed this regression and all four kinetic cases after the change.

The isolated serial Godunov iteration-2 diagnostic also passed after both
conservative admission fixes; the earlier endpoint-selection-only attempt
failed and remains in the evidence record.

The inventory change applies to all uses of the serial composite integrator,
including source ledgers and donor budget checks. Its last-bit rounding can
change a near-zero budget decision; no guard or tolerance was weakened. It is
not restricted to the admission option. The cell-update and reactor formulas
are unchanged; the diagnosed MOL1 input packets and local states were
bit-identical before this grouping fix.

All original exact-state and conservation assertions remain required. This
adversarial fixture reaches positive EOS-extrapolated temperatures near 17.7 K;
it is not physical combustion validation. No full-CFD, selected-mechanism matrix
or release qualification is claimed. The original no-capability kinetic mode
has not been fixed or qualified by these opt-in results.

Exact commands, build/test logs, source/binary hashes, independent reviews and
subsequent kinetic disposition belong to the phase3 integration evidence under
`docs/analysis/pelef_parallel_20261106/` in the GTT workspace.
