# Variable-grid reacting EB pilot, 2026-10-21

The 16-cubed root/fine-patch diagnostic completes two 50 ns root steps and
satisfies every unchanged field, species, element and numerical-energy-budget
limit against native held EOS. This establishes the finer-grid input path,
not spatial convergence or unchanged-PeleC qualification.

## Implementation and controls

The probe accepts a V2 input with root size 8, 16 or 32. The physical half-domain
fine box, ratio two, complete cell count, full-domain EB geometry, padded
geometry and covered-parent checks scale with that size. V1 retains its
8-cubed layout. The checker has an explicit `--mesh-study` mode that verifies
the native held-EOS executable and native activation-basis bundle, resizes
native grids, imports all initial cells and checks analytic fluid fractions.
Success is `diagnostic_passed`, with `qualification: false`; these captures
are refused by ordinary pinned-reference mode.

Both V1 and V2 8-cubed paths reproduce the earlier canonical-mechanism 25 ns
capture byte for byte, including all gas/UEINT/T/reaction-history fields and
budgets at zero and one step. Unsupported size, incomplete cell count and
diagnostic-capture reuse in pinned mode are refused. The runtime libraries
and prior executable/object files remain unchanged; new probe objects and
two linked executables live in an isolated directory.

## Completed finer-grid result

The 1 cm cube and planar wall at 3.3 mm are unchanged. The root and fine patch
each contain 4096 cells, with 4160 composite fluid cells (2624 coarse, 1536
fine). All 8192 initial gas/UEINT/T states match exactly. Both sides compute
their initial reaction histories independently for the chosen interval.

At 100 ns, the maximum temperature difference is `1.181633138e-5 K`, scaled
total-energy difference `3.048865040e-10`, and maximum z-momentum difference
`1.606224818e-7 kg/m2/s`. H2O mass-fraction difference is `2.430127933e-12`.
The relative hydrogen-consumption difference is `6.834390506e-7`;
maximum relative numerical-energy-budget residual is `5.587644133e-17`.
All three sampled times pass, with no tolerance changes.

Native held EOS took 12.05 seconds and the Fortran run 162.16 seconds.
These are diagnostic timings with runtime/FPE checks, not a release-build
performance comparison. GNU Fortran 13.3 and SUNDIALS 7.2.0 were reused.

## Evidence and next step

Artifacts: `build/eb-mesh-study-20261020/`; the completed report is
`runs/matched_eb_amr_6mo7eslw/checks.json`. `build.json` identifies the source,
compile/link commands and all linked objects/libraries. `compatibility.json`
records two successful compatibility runs and two input refusals;
`capture-rejection.json` records the third refusal. Source provenance and
regeneration, Python syntax, diff and project-contract checks pass.
`verification.json` records the completed identities and checks. There was no
full matrix, MPI/restart, install or physical-validation run.

The earlier 8-cubed native-activation trajectory used 100 ns root steps, so
it cannot isolate spatial differences from this 50 ns pilot. The next study
uses 50 ns on all three grids (8, 16 and 32), reusing this 16-cubed run.
Fine fields must be restricted with actual fluid volumes, excluding covered
coarse cells, before reporting spatial differences. Initial sampling error
must also be reported separately from the change during evolution.
