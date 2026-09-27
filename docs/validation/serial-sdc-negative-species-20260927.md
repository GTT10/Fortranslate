# Serial SDC accepted species — 2026-09-27

## Change and boundary

Implementation `38fa6776`, based on `091f0322` on `main`, prevents serial EB
tree SDC from publishing negative species densities that EOS trial-state
recovery otherwise accepts as roundoff. The common-step check runs after
reflux/redistribution, before any optional outputs or the candidate tree are
published. Recursive redistribution-transfer recovery also rejects negative
partial densities. Trial-state EOS tolerances remain unchanged; no species
clipping, normalization or golden-data changes are used.

This extends the [MPI acceptance repair](sdc-negative-species-20260927.md).
It does not connect serial MOL or fix the nested LF centroid positivity issue.

## Focused evidence

Two executables completed successfully, compiled with GNU Fortran 13.3.0,
`-g -fcheck=all -ffree-line-length-none`, against the rebuilt fixed-mechanism
`pelef_core` in `build/heptane-simple-mpi-debug` (Debug, tests disabled).
These are direct executable runs, not a CTest matrix or installed-binary audit.
The new test is registered as `unit_serial_sdc_species_3d` for tests-enabled
configurations; that CTest registration was not executed here.

1. `test_serial_sdc_species_3d`: a homogeneous periodic 2³ gas fixture receives
   a stage-1 auxiliary sink of `1e-8 kg/m³/s` in an initially absent species,
   with the same density source. A `1e-10 s` step would produce
   `-5e-19 kg/m³`. Both callbacks run. Before the fix, the negative-state
   rejection assertion fails. Afterward, the call fails transactionally and
   gas state, temperature, all five ledgers and minimum theta remain exactly
   unchanged. The sign-reversed source succeeds, producing `5e-19 kg/m³`
   and the matching auxiliary inventory. Two source cases pass.
2. Existing `test_eb_tree_full_3d`, default mode: common and recursive SDC,
   ancestor interpolation, rollback and explicit/implicit split cases pass.
   Reported SDC inventory residuals are `3.7650552248896492e-15` and
   `5.0897607707713024e-15`, respectively. Recursive negative-transfer
   injection was not separately exercised.

Early test drafts had an already-allocated callback buffer and an insufficient
geometry halo; these setup errors were corrected before the valid pre-fix
negative-state reproduction. The final before/after logs use the same test
source. Existing compiler warnings and executable-stack linker warnings remain.

Provenance preflight, project contract and `git diff --check` pass. No Python
source was changed, so Python syntax validation is not applicable.

## Reproduction artifacts

`build/serial-sdc-species/commands.txt` records exact build/run commands.
`before.log`, `after.log`, `tree.log`, and `build.log` retain the results.
Executable SHA-256 values:

| Artifact | SHA-256 |
| --- | --- |
| `before` | `c04b0e44bcecce052acbc6acde18c275273c773f322b84035ec48e8ded759ac1` |
| `after` | `55a0595390663070acc782cade783e731d330d0cb27f86cfa5adb76ea88bfdd4` |
| `tree_check` | `280d114ffd57138369efcf25fb3e47909b5a0ad5a4ba07755b9aec16766a36ee` |

The frozen running Jet binary was not rebuilt or stopped for this change.
