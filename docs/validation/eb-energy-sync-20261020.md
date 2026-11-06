# Reacting EB MOL energy synchronization, 2026-10-20

The ordinary reacting NoRedist MOL path now preserves the total/held-energy
difference left on a fine level by parent reflux until that level advances.
The serial and MPI drivers validate the actual held-energy EOS at clock
entry; valid states no longer incur an extra total-energy overwrite there.
Paired synchronization resets only the closing level. Other levels undergo
held-EOS validation. Invalid held states still enter the existing dual-energy
repair and correction ledger; an unsuccessful repair refuses the transaction.
Predictor/corrector repairs remain in place. CTU, nonreacting and prepared
StateRedist paths retain their prior reset policy.

Pressure recovery continues to use actual stored UEINT, including under
large kinetic/rotating-energy cancellation. This change does not replace
held-energy pressure with total-minus-kinetic pressure or change tolerances.
Ordinary native replay entry follows the same admission rule. Common-clock
closure and recursive MPI closure after recovery use the same level scope.

## Diagnosis

The preceding [transport correction](eb-transport-centroid-20261020.md)
left a 300 ns temperature difference of `1.966751177e-3 K` and 19 failing
field limits. Native `PeleC.cpp` applies reflux to coarse and fine states,
then `post_timestep` computes temperature on its own level. At 100 ns the
saved native fine state retains a `0.167289543 J/m3` gap between total-minus-
kinetic energy and stored internal energy. The previous Fortran closure
reset all levels and its next clock entry repeated that reset.

Three fresh diagnostic interventions used copied modules, with no alteration
of pinned sources, the saved native reference, or acceptance limits:

| Intervention | 300 ns temperature error (K) | Total-energy error (J/m3) | Failing limits |
| --- | ---: | ---: | ---: |
| Reset only closing level, retain old entry resets | 1.957469459e-3 | 0.300556724 | 19 |
| Also omit the three old-clock entry resets | 3.848490223e-5 | 0.02020013539 | 12 |
| Additionally use total-derived hydro pressure | 2.908427496e-5 | 0.01183452041 | 12 |

Only the timing change is implemented. The third diagnostic still fails
the reference gate and would change the held-energy cancellation contract.
The production implementation validates entry states and retains repair
fallback rather than unconditionally skipping admission.

The checker now reports `rho_E - rho_e - |momentum|^2/(2*rho)` in SI units
for each composite fluid level, excluding covered coarse cells. The report
includes the maximum, signed worst-cell value and volume integral. This is
diagnostic metadata, explicitly **not an acceptance gate**. A small analytic
check covers both gap signs, volume weights, empty input, invalid density and
the captured 352 coarse plus 192 fine composite cells.

## Reproduction and evidence boundary

Artifacts: `build/eb-energy-sync-20261020/`. Base: dirty `main` at
`0a0e00ebe9138145fbb0e0e745665d0d768df009`. `before/`, `increment.patch`,
source copies and hash manifests separate this increment from prior work.
The native capture is reused from
`build/tgreact-eb-amr-reference-20261019/runs/matched_eb_amr_62p3fjsx/reference`;
native code was not rerun. The case uses three 100 ns root steps, ratio-two
subcycling, MOL2, matched Simple LiDryer transport and CVODE, with the same
`eta2=2`, negative-energy policy and field/inventory thresholds on both sides.

Build: `build/reacting-hot-diagnosis-20261019/registration`, GNU Fortran
13.3.0 Debug with bounds/runtime and invalid/zero/overflow checks, Open MPI
4.1.6, CMake 3.28.3 and SUNDIALS 7.2.0. No build configuration changed.
MPI uses `/usr/bin/mpiexec`, `OMP_NUM_THREADS=1`, `HWLOC_COMPONENTS=-gl,-x86`.

Diagnostic runs are `mode-1/matched_eb_amr_l2170h3i`,
`mode-2/matched_eb_amr_3vt2be7m`, and
`total-primitive/matched_eb_amr_pn5r3mtn`. Exact commands and copied-module
identities are retained. The first two checker copies predate only the
top-level diagnostic-definition metadata; `checker-at-diagnostics.py` matches
their recorded SHA. Two premature launches failed before the diagnostic
binary existed and are recorded in `diagnostic-startup-failures.json`.
The first production build failed on added lines exceeding Fortran's line
length limit; `build-production-line-length-failure.log` retains the failure.
These startup/build failures are not numerical results.

Full configuration-matrix, installed-binary and physical/ignition validation
are outside this focused change. Pinned nonuniform EB-AMR trajectory parity
remains open.

## Completed production verification

The final production build completed successfully (`build-production-final.log`).
Execution session `79476` completed all three checks; results were inspected
on continuation. Exact build/check commands, source manifests and executable
hashes are retained in `verification.json` and `check-commands.json`.

- `unit_eb_mol_hydro_3d` passed in 2.86 s, including actual held-energy and
  rotating cancellation coverage.
- `unit_mpi_eb_native_mol_partial_recovery_3d_3` passed in 12.72 s.
- Hot CVODE subcycling, reaction-off control and changed-rank restart passed
  in 68.47 s: five successful runs and two expected invalid-input refusals.
  Intermediate/final gas exports match byte-for-byte across serial, MPI and
  2-to-3-rank restart. Evidence:
  `application-runs/reacting_restart_yov6pz1z/checks.json`.
- Mechanism provenance with regeneration, Python syntax, diff whitespace,
  project contract and the analytic energy-diagnostic check passed. No skips.

The production fixed-reference probe completed all three root steps but
**fails 12 unchanged field limits**, down from 19 before this change.
All first-step field limits now pass. At steps two and three, density,
three momentum components, total energy and held energy still fail.
Temperature and species mass fractions pass at all compared times. Evidence:
`reference-runs/matched_eb_amr_pnnsdcao/checks.json` (comparison exit 1;
the Fortran solver itself exits 0).

| Maximum absolute difference at 300 ns | Before | Production after |
| --- | ---: | ---: |
| Density (kg/m3) | 3.702090356e-8 | 5.972077999e-9 |
| x-momentum (kg/m2/s) | 3.201288544e-5 | 4.764517059e-6 |
| z-momentum (kg/m2/s) | 1.160216988e-5 | 2.197937003e-6 |
| Total energy (J/m3) | 0.457020083 | 0.02020013545 |
| Held energy (J/m3) | 0.300697281 | 0.02018125047 |
| Temperature (K) | 1.966751177e-3 | 1.098624132e-4 |

Production temperature differs from the copied-module diagnostic because
successful held-EOS validation refreshes temperature even when it preserves
gas and UEINT. At 100 ns, every exported gas/UEINT value matches diagnostic
mode 2 exactly while temperature differs by up to `8.067916315e-5 K`.
At 300 ns, their energy values differ by at most `1.75e-10 J/m3` and
temperature by `7.46e-5 K`. The smaller diagnostic temperature error must
not be reported as the production result. `production-comparison.json`
records these comparisons.

The production fine-level total/held gap at 300 ns is `0.156514949 J/m3`,
versus native `0.156528647 J/m3`; the coarse gap remains at roundoff.
Maximum Fortran mass/element inventory drift is `5.98e-14`; maximum energy
budget residual, including the recorded numerical correction, is `9.38e-17`.
Hydrogen-consumption relative error is `2.24e-5`. These are conservation and
short-time chemistry evidence, not complete trajectory parity. The remaining
field maxima lie on the fine level. Further work must separate its second-
clock spatial/forcing discrepancies from the retained pressure-policy
difference, without changing held-energy cancellation safety or thresholds.

## Second-clock isolation

A copied probe omitting only startup reaction-history initialization imports
the saved 100 ns gas/T/UEINT and reaction history, then advances one 100 ns
root interval with the final production libraries. Absolute output time is
offset by 100 ns for this diagnostic comparison; the case has no time-varying
boundary or external forcing. The reference capture is unchanged.
Artifacts: `build/eb-second-clock-20261020/`, including the copied probe,
`run_diagnostic.py`, exact compile/link/run commands and hashes.

| Imported 100 ns state / history | 200 ns x-momentum error (kg/m2/s) | Total-energy error (J/m3) |
| --- | ---: | ---: |
| Fortran / Fortran | 3.335217856e-6 | 0.01331079935 |
| Native / native | 3.299236876e-6 | 0.01326166201 |
| Fortran / native | 3.335218110e-6 | 0.01331078660 |

All three diagnostic runs complete in 17.1–17.5 s. The Fortran/Fortran control
reproduces the uninterrupted production endpoint exactly in all compared
gas/T/UEINT fields. Replacing the reaction history alone has negligible
effect. Importing both the native state and its history also leaves nearly
the same fine-level discrepancy. Thus accumulated startup/history error is
not the principal cause; the discrepancy is generated within the next
interval. This narrows subsequent work to that interval's spatial/forcing
evaluation and retained energy/pressure policy. It does not qualify native
restart or establish complete parity.

The first diagnostic compile omitted the AMR module include directory and
failed before running. The corrected compile uses the production target's
include paths; `build-missing-include.log` retains the failure. No production
numerical source or acceptance threshold changed during this isolation.
