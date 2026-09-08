# Current status: 0.245.0 integration candidate

This page separates present capability, dated reported results, and independent
integration checks. It supersedes global-status interpretations of older
milestone prose, not the original test evidence or implementation details.

## Supported boundary

| Area | Implemented subset | Still open |
| --- | --- | --- |
| Reactive flow and transport | 1D/2D paths and periodic regular 3D NASA7 mixtures, cell-local chemistry, molecular transport | General production 3D boundary/coupling qualification |
| 3D AMR | Static, strictly interior, periodic, ratio-two two-level hierarchy; fixed/selected chemistry/transport; sparse MPI; transport-aware restart | Dynamic topology, multiple general patches/levels, physical coarse boundaries |
| 3D embedded boundary | Bounded planar, single-level hydro/chemistry/transport/restart paths | General geometry and integrated 3D EB-AMR-MPI qualification |
| Chemistry | Fixed H2/O2 and configure-time selected supported bundles; optional official-Fortran CVODE for selected 0D | General detailed-fuel ingestion/validation and stiff CFD integration |
| Parallelism/output | Rank-parity and changed-rank restart for qualified subsets | Production scaling and scalable/distributed I/O |
| LES and spray | Not claimed complete | Validated LES, particle transport/evaporation/coupling and breakup |

Passing conservation, rollback, byte-parity, or short-time H2/O2 tests does not
establish physical ignition/spray validation or complete PeleC equivalence.
License selection remains an owner decision, not part of this repair.

## Evidence ledger

The checkpoint `509d07929d52d9e570a6af1f6ec7e5836bd4fe5d` saved the 0.245.0
worktree. Its initial PR #102 CI failed at Configure. The independently
reproduced cause was two missing spaces in the vendored YAML description;
restoring the exact Cantera 3.2.0 package source leaves the chemistry bundle
and generated Fortran unchanged. See [the repair record](validation/checkpoint-integrity-20260908.md).

**Reported historical local qualification:** the
[2026-09-08 handoff](handoff-pelef-0.245.0-20260908.md) records 4570/4570 across
Debug 524, Release 524, Debug+Cantera 529, MPI Debug 694, MPI Release 694,
CVODE Debug 533, CVODE Debug+Cantera 539, and CVODE Release 533. It also reports
a tests-disabled 42/42 ELF install audit, separate tests-enabled 23/23 installed
restart tests, and installed application/checker/byte comparisons of 32/32,
8/8, and 12/12. These are reported results tied to its frozen worktree hash and
local evidence paths; this repair did not retrieve those raw server directories
or independently rerun that eight-configuration matrix. Do not include the
interrupted 227-test MPI run or the excluded literal-prefix exploratory output.

**Independent integrity qualification:** commit
`566010ffc77a6de25d452d084ff61fa2bb533d95` passed
[mechanism-integrity run 34208154353](https://github.com/GTT10/Fortranslate/actions/runs/34208154353),
including 14 checker regressions, exact source hashing, and full Cantera bundle
regeneration. This proves the source/bundle repair, not the complete numerical
release matrix. Final serial/MPI/optional-backend integration status belongs to
[PR #102](https://github.com/GTT10/Fortranslate/pull/102) and dated follow-up
validation records, with their tested commit and configuration.

The earlier focused 10/10 serial and 23/23 MPI results remain in
[validation/0.245.0.md](validation/0.245.0.md); the committed older install report
also retains its own date and source provenance. Neither should be silently
overwritten with the later handoff counts.

## Immediate work

Finish repaired-checkout serial/MPI/install verification and preserve the
remaining optional-backend evidence before closing integration. Then follow
[the roadmap](completion_roadmap.md), one acceptance-tested increment at a time.
