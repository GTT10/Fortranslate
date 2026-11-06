# Covered coarse energy after reacting MOL closure, 2026-10-20

The preceding [energy synchronization change](eb-energy-sync-20261020.md)
left twelve reference field limits open. A second-clock state/history
intervention showed that those discrepancies are generated within the next
interval, even when that interval starts from native gas and reaction history.

## Stage diagnosis

Fresh copied-module captures record all eighteen spatial evaluations in two
root MOL2 steps, including full available halos, gas, UEINT, pressure and
separate Fortran hydro/transport rates. The native capture changes only
diagnostic output in a copy of `Diffusion.cpp`; pinned source files are intact.
All three native and all three Fortran time-level CSV exports are byte-identical
to the preceding production captures. Thus instrumentation does not change
the measured trajectory. Evidence: `build/eb-second-clock-spatial-20261020/`,
especially `comparison.json`, `capture-verification.json`, source copies,
hashes and exact build/run commands. The inherited Fortran executable name
`probe_total_primitive` is historical: this capture uses the current held-EOS
hydro source without any pressure substitution.

At the second root entry (spatial evaluation 9), the largest coarse total-
energy difference is `0.282184763 J/m3`, at covered coarse cell `(4,4,4)`.
It becomes a `0.141150546 J/m3` total/held-energy difference in that cell at
the second corrector (evaluation 11), changing pressure by `0.0432053 Pa`
and the subsequent fine halo. In the previously exported Fortran state, the
covered cell's `E - H - kinetic` gap is `0.281960749 J/m3` at 100 ns and
`0.315052245 J/m3` at 200 ns. These covered cells do not enter composite
inventory checks but do supply interpolation and coarse evolution.

Pinned `PeleC::post_timestep` calls `avgDown()` and then `computeTemp()` on
the closing level; `computeTemp` first resets energy, including covered coarse
cells. PeleF restricted the gas/held families after its leaf reset, leaving
covered total energy inconsistent with held energy. Redundant outer-driver
restriction could also undo a corrected covered state before its next use.
This defect is separate from the documented total-derived versus held-energy
pressure policy on fine leaves.

## Implementation

Ordinary reacting NoRedist MOL closure now resets covered fluid cells on the
closing level after restriction, in serial and MPI. The existing dual-energy
policy and EOS perform that reset; rotating cancellation handling is unchanged.
The operation runs on the private candidate before publication. Composite
energy and held-energy ledgers exclude covered cells, and the closure checks
that this second reset contributes zero composite energy change.

The common and recursive driver entry/exit paths no longer repeat restriction
for this mode after a clock has already performed its paired closure. Other
integrators and redistribution modes keep their existing behavior. Fine leaves
retain the post-parent-reflux total/held gap from the preceding timing fix.

The fixed-case production probe now independently checks the total/held/kinetic
identity in covered coarse fluid cells after each root step, using its fixed
`eta2=2` policy. This check would detect the preceding `0.282 J/m3` gap even
while the full reference trajectory gate remains open. A focused artifact
check also checks that fine-level gaps survive and the composite energy
budget remains closed; no reference tolerances are changed.

Artifacts: `build/eb-covered-energy-20261020/`. Base: dirty `main` at
`0a0e00ebe9138145fbb0e0e745665d0d768df009`, preserving the preceding numerical
changes. `before/`, `after/`, `increment.patch` and `source-identity.json`
identify this increment. Build configuration remains GNU Fortran 13.3.0 Debug,
runtime/bounds/FPE checks, Open MPI 4.1.6 and SUNDIALS 7.2.0 in
`build/reacting-hot-diagnosis-20261019/registration`.

## Completed verification

The final build and execution session `37667` completed successfully, except
for the unchanged fixed-reference field gate described below. Mechanism
provenance with regeneration passed. Exact commands, source and executable
hashes, and final static checks are retained in `verification.json`.

- `unit_eb_mol_hydro_3d` passed in 2.86 s.
- `unit_mpi_eb_native_mol_partial_recovery_3d_3` passed in 12.41 s.
- Hot reacting subcycling, the reaction-off control and changed-rank restart
  passed in 67.68 s: five successful runs and two expected refusals. Gas
  snapshots match byte-for-byte across serial/MPI and 2-to-3-rank restart.
  Evidence: `application-runs/reacting_restart_pbrxonni/checks.json`.
- The production probe's covered-energy assertion passes after each root
  step. The 32 covered coarse fluid cells have a maximum energy gap of
  `1.45e-11 J/m3`. Fine-level gaps remain `0.1673`, `0.1616` and `0.1565 J/m3`,
  confirming that the fine leaves were not reset by this correction.
  Evidence: `covered-energy.json`.

The three-step reference probe completes, but **twelve original field limits
still fail**. At 300 ns, the density error is `5.968e-9 kg/m3`, x-momentum
`4.765e-6 kg/m2/s`, total energy `0.0202096 J/m3`, held energy
`0.0201809 J/m3` and temperature `1.09801e-4 K`. These are nearly unchanged
from the preceding production result. The covered-parent inconsistency is
therefore corrected, but is **not the principal cause** of the remaining
trajectory discrepancy. The maximum corrected composite energy-budget
residual is `1.19e-16`; hydrogen-consumption relative error is `2.24e-5`.
Evidence: `reference-runs/matched_eb_amr_z41j_1e2/checks.json` (comparison
exit 1, Fortran solver exit 0). No acceptance thresholds were changed.

The stage diagnostic's first transport-module compile lacked the EB include
directory; `build-missing-include.log` retains that failure. After using the
production include directories, both captures and exact-export checks pass.
No full configuration matrix, installed-binary audit or physical validation
has been run for this increment.
