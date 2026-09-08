# PeleF repository instructions

## Start here

Read `README.md`, `docs/current_status.md`, and the relevant case/module.
Use `docs/README.md` as an index; do not load every historical milestone into
context. Inspect `git status`, the branch, and the intended base before editing.
`docs/handoff-pelef-0.245.0-20260908.md` is a dated historical handoff, not a
standing instruction to reuse its old worktree or overwrite current user work.

## Change boundaries

- Keep one purpose per change. Separate numerical changes from formatting,
  build-system refactoring, and qualification-record edits.
- Preserve user changes and running simulations. Do not reset, clean, force
  push, delete branches, or stop unrelated jobs without explicit authorization.
- Use current `main` for new work. PR #102 is integrated and historical PR #11
  is closed as superseded; do not reopen their accidental PR stack.
- Spray/SGS support is opt-in and experimental. Read `docs/spray.md` before
  extending it. Do not label periodic small-grid conservation tests as nozzle,
  diesel ignition, AMR-particle or validated LES results.
- Read generated-file headers. Change a mechanism's source and regenerate with
  the pinned importer; never hand-edit generated rates or bypass SHA checks.
  `mechanisms/h2o2_cantera.yaml` is byte-pinned upstream data: preserve its
  intentional whitespace, including two spaces on line 4.
- Preserve transactional failure/rollback, species and element conservation,
  rank-neutral restart, and the fixed/selected mechanism boundary. Do not
  loosen tolerances or replace golden hashes merely to turn CI green.

## Verification and reporting

Run the provenance preflight, relevant unit/application tests, Python syntax,
`git diff --check`, and `tools/check_project_contract.py` for touched areas.
Use clean build directories for compiler/MPI/optional-backend changes.
Record exact source identity, commands, toolchain, test counts, failures,
skips, and artifact locations. Report only completed checks as passed.

A focused test, full configuration matrix, installed-binary audit, external
solver parity, and physical validation are different evidence classes.
Keep reported historical results separate from independently rerun results.
Do not mix tests-enabled/disabled install prefixes or count interrupted runs.

## Documentation ownership

`README.md` is the entry point; `docs/current_status.md` is the present scope;
`docs/completion_roadmap.md` and Issues own future work; design decisions explain
why; `docs/validation/` owns dated evidence. Existing long milestone inventories
are historical detail. Do not duplicate each implementation log across all
of these files or increment the numerical milestone for housekeeping alone.
