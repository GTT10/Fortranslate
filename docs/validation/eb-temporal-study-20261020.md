# Reacting EB temporal diagnostic, 2026-10-20

The bounded time study completes all three schedules with conserved mass and
elements and a closed PeleF numerical-energy budget. The native held-EOS
diagnostic agrees with PeleF at every sampled time through 300 ns for 50 ns
and 25 ns root steps. At 100 ns root steps, O2 and H2O exceed the unchanged
mass-fraction limit starting at 1.0 microsecond. This is an unresolved longer
trajectory difference, not a passing qualification.

## Scope and identity

`tools/check_tgreact_eb_temporal_study.py` reuses the completed production
capture and separately identified [native held-EOS diagnostic](eb-held-policy-20261020.md).
Production numerical sources, binaries, case geometry, initial gas/UEINT/T,
Simple transport fits, mechanism and tolerances remain unchanged. Initial
reaction histories are computed independently for each schedule's time step.
The native diagnostic still differs from the pinned spatial EOS algorithm;
the original twelve pinned-reference failures remain open.

The mesh remains an 8-cubed root plus an 8-cubed ratio-two fine patch, with
544 composite fluid cells. This study does not establish spatial convergence,
ignition, experimental agreement or MPI behavior at these longer times.

| Root step | Endpoint | Steps per solver | Fortran runtime | Native runtime | Diagnostic field exceedances |
| --- | --- | ---: | ---: | ---: | ---: |
| 50 ns | 300 ns | 6 | 97.71 s | 5.74 s | 0 |
| 25 ns | 300 ns | 12 | 182.72 s | 7.28 s | 0 |
| 100 ns | 1.2 microseconds | 12 | 214.71 s | 24.44 s | 6 |

Times are wall times with three concurrent schedules and one thread per
executable, not a controlled performance comparison. The existing 100 ns,
300 ns production and native captures supply the coarsest time resolution;
they were not rerun. The three new schedules contain 30 root steps per solver
and 33 snapshots per solver, including their initial states.

## Observed temporal differences

At 300 ns, use composite fluid-volume-weighted RMS differences between
100/50 ns and 50/25 ns solutions. The observed order is log2 of their ratio.
An independent calculation applies the same norm to native held-EOS snapshots.

| Field | PeleF observed order | Native held-EOS observed order |
| --- | ---: | ---: |
| Density | 1.125488 | 1.125488 |
| x-momentum | 1.188849 | 1.188849 |
| Total energy | 1.153778 | 1.153777 |
| Held energy | 1.152922 | 1.152922 |
| Stored temperature | 1.207223 | 1.207196 |

Differences decrease, but this interval does **not** demonstrate second-order
time accuracy. The same trend in both implementations means it is not evidence
of a Fortran-only time-integrator defect. Startup, AMR coupling and asymptotic
behavior remain unresolved. The table is an observed three-resolution trend,
not a demonstrated asymptotic order.

## Longer trajectory

At 1.2 microseconds, maximum O2 and H2O mass-fraction differences are
`5.419574964e-7` and `4.469964739e-7`, respectively, versus the unchanged
`1e-7` limit. Each exceeds that limit at steps 10, 11 and 12. All other field
limits pass through this endpoint. Temperature differs by `2.299841849e-4 K`;
scaled total-energy difference is `3.420620610e-9`.

PeleF consumes 0.7226098% of the initial hydrogen; native held EOS consumes
0.7226810%. Their relative consumption difference is `9.84629155e-5`, within
the `1e-3` limit. Across all schedules, maximum mass/element drift is
`5.9730e-14` for PeleF and `4.3299e-15` for native; maximum relative PeleF
energy-budget residual is `1.2912e-16`. Species remain nonnegative and their
sum closes. These checks do not override the six field failures.

The next controlled investigation is the chemical source mismatch. Source
inspection finds that pinned CEPTR uses `Rc=1.98721558317399615845 cal/mol/K`
for activation temperatures, whereas the canonical imported Cantera mechanism
uses the SI calorie conversion. All 17 nonzero native activation coefficients
agree with the former conversion to `2.11e-15` relative. That source-level
difference is a hypothesis for the accumulated chemical discrepancy; a
separately generated compatibility variant is required to test its effect.

## Reproduction and checks

Artifacts: `build/eb-temporal-study-20261020/`, with complete schedules under
`runs/temporal_pqfpp6bu/`. `checks.json` preserves field failures, inventories,
consumption, command times, initial/source/binary identities and convergence
norms. `native-temporal-comparison.json` records independent native norms and
input snapshot hashes. `resume.json` contains the exact driver command;
`verification.json` records completed checks and source identity.

The reused toolchain is GNU Fortran 13.3, Debug with runtime/FPE checks,
OpenMPI 4.1.6 and SUNDIALS 7.2.0. LiDryer provenance and regeneration passed
before the runs. Python syntax, project contract and `git diff --check` pass.
No unit suite, full configuration matrix, installation audit or physical
validation was rerun. All six solver executions completed; none was counted
from interrupted output. Diagnostic exceedances remain explicitly reported.
