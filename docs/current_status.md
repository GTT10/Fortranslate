# Current status: PeleF 0.246.0

This page separates implemented capability, dated reported results, and
independent checks. The preceding integration is
recorded in [PR #102](https://github.com/GTT10/Fortranslate/pull/102). The new
opt-in spray scope and qualification boundary are in
[the spray guide](spray.md) and [0.246.0 validation](validation/0.246.0.md).
Older milestone prose remains historical detail, not a current global status.

## Supported boundary

| Area | Implemented subset | Still open |
| --- | --- | --- |
| Reactive flow and transport | 1D/2D paths and periodic regular 3D NASA7 mixtures, cell-local chemistry, molecular transport | General production 3D boundary/coupling qualification |
| 3D AMR | Static, strictly interior, periodic, ratio-two two-level hierarchy; fixed/selected chemistry/transport; sparse MPI; transport-aware restart | Dynamic topology, multiple general patches/levels, physical coarse boundaries |
| 3D embedded boundary | Bounded planar, single-level hydro/chemistry/transport/restart paths | General geometry and integrated 3D EB-AMR-MPI qualification |
| Chemistry | Fixed H2/O2 and configure-time selected supported bundles; official-Fortran CVODE for selected 0D and opt-in regular 3D spray cells; FFCM1_Red methanol fixture | 2--32 species; production diesel-fuel mechanisms and sparse/performance qualification |
| Parallelism/output | Rank-parity and changed-rank restart for qualified subsets | Production scaling and scalable/distributed I/O |
| LES and spray | Opt-in periodic regular 3D single-component parcels; drag/heating/evaporation, conservative two-way coupling, deterministic cone/scheduled injection and mass-CDF size distributions, deviatoric Smagorinsky, coupled restart | Multicomponent/boiling/breakup/collision, particle MPI/AMR/EB, physical boundaries, calibrated properties and experimental validation |

Conservation, rollback, byte-parity, and short-time H2/O2 tests do not establish
physical ignition/spray validation or complete PeleC equivalence. Choosing a
license remains an owner decision.

## Repository consolidation implemented

The README is an entry point rather than a chronological implementation log.
Contributor instructions, AGENTS.md, the documentation index and PR template
separate current capability, historical evidence and future work.

The 2613-line root CMake and 9138-line test registration were split into five
root target modules and nineteen test modules, preserving include order and
directory scope. Generated target graphs, install rules and CTest commands
matched in seven clean configurations; see [the layout record](validation/build-layout-20260909.md)
and [CMake navigation](../cmake/README.md).

Hosted workflows cover source provenance/regeneration, complete serial and MPI
suites, selected 0D CVODE with pinned SUNDIALS 7.2.0 and Cantera 3.2.0, and a
separate tests-disabled Release installation. Source identity, dependency
versions, JUnit and logs are retained. An implemented workflow is not itself a
passed qualification; completed results must be associated with its tested SHA.

Historical draft #11 was reviewed and closed as superseded: its 17 numerical,
mechanism, application/case and verification paths were byte-identical at the
later integrated full-chemistry commit. Its unique decision note is preserved
[in the historical record](history/pr11-implicit-chemistry-decision.md).
No historical branch was deleted. The audited ancestor-only bulk cleanup was
blocked by the connection's safety restriction, not performed silently.

## Evidence ledger

The saved checkpoint is `509d07929d52d9e570a6af1f6ec7e5836bd4fe5d`.
Initial PR CI failed because its vendored YAML description lacked two spaces.
Restoring the exact Cantera package bytes repaired source provenance without
changing the mechanism bundle or generated Fortran; see [the integrity repair](validation/checkpoint-integrity-20260908.md).

A later targeted numerical review found a false-success return in NASA7 molar
conversion overflow guards. Commit `bb2980e52c2e6a4bfed353dd280fe1f0e21f4893`
corrects the failure flag without changing physical formulas or tolerances.
New cp/enthalpy/entropy rejection cases fail on the parent implementation and
pass after the fix; seven focused GNU14 Debug tests passed. See [red/green evidence](validation/nasa7-conversion-status-20260909.md).
This correction is intentionally separate from the behavior-preserving layout.

**Reported historical local qualification:** the [2026-09-08 handoff](handoff-pelef-0.245.0-20260908.md)
reports 4570/4570 across Debug 524, Release 524, Debug+Cantera 529, MPI Debug 694,
MPI Release 694, CVODE Debug 533, CVODE Debug+Cantera 539, and CVODE Release 533.
It also reports tests-disabled 42/42 ELF, separate tests-enabled 23/23 installed
restart checks, and installed application/checker/byte comparisons of 32/32,
8/8 and 12/12. These belong to that frozen worktree and its local evidence;
the raw server directories were not retrieved here. Interrupted 227-test and
excluded literal-prefix exploratory results are not successful qualification.

**Independent source-integrity qualification:** commit
`566010ffc77a6de25d452d084ff61fa2bb533d95` passed
[run 34208154353](https://github.com/GTT10/Fortranslate/actions/runs/34208154353),
including 14 checker regressions, exact hashing and full bundle regeneration.
Final serial/MPI/CVODE/clean-install results belong to PR #102 and their own
source identities. They must not be replaced by older source successes.

The earlier focused 10/10 serial and 23/23 MPI results remain in
[validation/0.245.0.md](validation/0.245.0.md). Historical evidence is preserved,
not overwritten or added together as though all configurations were identical.

## Remaining work

[The roadmap](completion_roadmap.md) distinguishes further solver development
from repository consolidation. General 3D AMR/EB, production detailed-fuel stiff CFD,
full LES/spray qualification and external physical validation remain genuine engineering work,
not features implicitly completed by housekeeping or large test counts.
