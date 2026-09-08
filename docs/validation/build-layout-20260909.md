# CMake layout and hosted qualification, 2026-09-09 (JST)

## Source and scope

Before: commit `bf74884043a8cf25aece292604dd5d88035831bf`, tree
`795f7c7021ddf650008506f02545babfa27a6d22`.
After the layout-only change: commit
`e582c5dee11ac50f5cdb965d8e022287ccc0bd91`, tree
`8577536b78a62975163ecb4f0bcf3607088f5381`.

The 2613-line root and 9138-line test CMake files were separated into five
root target modules and nineteen test modules. The expanded contents are
byte-identical to the original files; include order and directory scope are
preserved. No Fortran, input case, mechanism, golden result, or tolerance was
changed. The one-time migration tools live only on the isolated maintenance
branch and are not added to the integration tree.

## Completed generated-contract comparisons

[Workflow 34252027782](https://github.com/GTT10/Fortranslate/actions/runs/34252027782)
completed successfully. Every configuration below was configured both before
and after the change using a clean build directory. All generated CTest
commands/properties, target dependency graphs, and install rules match after
normalizing source/build/prefix paths and backtrace-only CMake metadata.

| Configuration | Registered tests before / after |
| --- | --- |
| Debug, fixed | 524 / 524 |
| Release, fixed | 524 / 524 |
| Debug with Cantera | 529 / 529 |
| MPI Debug with selected H2/O2 | 706 / 706 |
| MPI Release with selected H2/O2 | 706 / 706 |
| Selected serial Release, tests disabled | 0 / 0 |
| Selected MPI Release, tests disabled | 0 / 0 |

These numbers are registration counts, not a claim that this comparison ran
all numerical tests. The 706-test selected MPI configuration is distinct from
the historical 694-test MPI configurations. The provenance unit tests and
project-contract checker also passed in this job.

[Artifact 10066398999](https://github.com/GTT10/Fortranslate/actions/runs/34252027782/artifacts/10066398999)
contains source commit/tree, the result source archive, configuration logs and
`layout-verification.json`. Archive SHA-256:
`f2aa8c3a5bd63c5f2333805a59ac9dcb5920ec9fd2b258026926a7b4e61c5d64`.
The archive was independently downloaded and its SHA-256 verified.

## Hosted numerical qualification added after the layout split

The integration follow-up adds an optional-backend/clean-install workflow:
SUNDIALS 7.2.0 is checked out at immutable commit
`71a4cc9ad5e7bc8b4e33a1ca9795b4e96883f9a6`, compiled with its official static
Fortran interface and double precision, and used by Debug/Release selected
0D CVODE builds with Cantera 3.2.0 parity enabled. Each configuration runs its
complete registered suite and preserves JUnit, logs and dependency identity.

A separate clean runner builds MPI Release with tests, CVODE and Cantera
reference tests disabled, checks the existing 42-program ELF installation
contract, and runs the existing installed application/restart smoke checks.
It does not reuse a tests-enabled install prefix.

Serial CI now retains source snapshots and JUnit/log artifacts. MPI CI retains
all previous commands and numerical conditions; its overall job budget is
raised from 35 to 90 minutes, and evidence retention from 14 to 30 days.
The per-test timeouts and numerical acceptance criteria are not relaxed.

Results for these new workflow jobs must be read from the current tested SHA
in [PR #102](https://github.com/GTT10/Fortranslate/pull/102). Adding this workflow
is not itself a successful numerical qualification. Historical 4570-test
worktree evidence remains a separate record.
