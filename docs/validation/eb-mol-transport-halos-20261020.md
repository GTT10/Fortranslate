# Held EB MOL transport halos, 2026-10-20

Ordinary held MOL now evaluates molecular transport with the same prepared
six-layer state/geometry neighborhood used by its hydro operator. Connected
coarse-to-fine interpolation is enabled for both operators. Serial and MPI
use the existing prepared transport kernels, without applying StateRedist
when the caller selects NoRedist. The CTU and gas-only default paths retain
their existing behavior.

## Diagnosis

Artifacts are under `build/eb-mol-primitive-energy-20261020/`.
The preceding full-transport comparison had a 300 ns temperature difference
of `2.781487e-3 K` and x-momentum difference of `6.557755e-4 kg/m2/s`.

A separately compiled copy of `eb_mol_hydro_3d_mod.F90` changes only the
spatial primitive and wall-pressure energy input to total-minus-kinetic.
The normal source and executable are unchanged by this diagnostic. Its
three-step comparison completes but still fails 20 field limits: the
maximum temperature difference is `2.781491e-3 K`, effectively unchanged.
Evidence: `diagnostic-runs/matched_eb_amr_n7ttu5at/checks.json`.

A second diagnostic compares initial spatial rates and records each MOL stage
in a one-step transport-off run. At the first fine predictor the maximum
x-momentum input difference is `4.02e-9 kg/m2/s`; at the second substep's
predictor it is `6.62e-8`. This distinguishes the small hydro-stage difference
from the much larger full-transport endpoint difference. It does not claim
transport-off endpoint qualification. The initial total-energy primitive
intervention reduces instantaneous rate differences, but that observation
alone does not improve the coupled trajectory.

Inspection found ordinary transport still using the face-donor halo path,
while held MOL hydro used the complete, spatially interpolated neighborhood.
The prepared transport implementation already supports full same-resolution
stencils. Its face-centroid transport correction remains missing (see below). The new connection uses it for ordinary held
MOL and opts the StateRedist MOL path into the identical coarse interpolation
as hydro. The existing exact old-halo equality check remains in place.

The independent held energy, its EOS admission, and its failure/rollback
behavior are preserved. In particular, no total-energy substitution is made
in production: the existing rotating high-cancellation pressure check remains
applicable. Molecular transport and wall exchange retain their conservative
face accounting. The new MPI interpolation control is collectively agreed
before any donor exchange, and gas fields remain owner-local.

## Reproduction

The working base is dirty `main` at
`0a0e00ebe9138145fbb0e0e745665d0d768df009`. `before/`, `source/`,
`increment.patch`, `before-identity.json`, and `source-identity.json` retain
this increment separately from existing changes. Diagnostic compiler/link
commands, original source copies, stage dumps, and comparisons are retained
in the same artifact directory. The pinned native source, binary, case hashes
and numerical limits are unchanged.

Build: `build/reacting-hot-diagnosis-20261019/registration`, GNU Fortran 13.3.0
Debug with runtime/floating-point checks, Open MPI 4.1.6 and CMake 3.28.3.
`check-commands.json` and `check-results.json` retain the focused MPI and
application commands and outcomes. The fixed-reference run retains its own
commands, executable hashes and outcomes in `checks.json`. MPI uses
`/usr/bin/mpiexec`; `OMP_NUM_THREADS=1`.

## Completed checks

- Mechanism provenance preflight with regeneration passed:
  `build/cantera-venv/bin/python tools/check_mechanism_provenance.py --bundle mechanisms/lidryer.json --regenerate`.
- Incremental build of the reference, held kernel, MPI partial recovery and
  selected serial/MPI applications passed; `build.log` and `final-build.log`.
- Two focused CTests passed, with no skips: `unit_eb_mol_hydro_3d` (2.94 s)
  and `unit_mpi_eb_native_mol_partial_recovery_3d_3` (12.72 s).
  The held kernel includes the independent rotating high-cancellation pressure
  and late held-EOS refusal/rollback cases.
- The hot CVODE subcycling application passed in 59.06 s: five successful
  executions (serial, reaction-off, MPI2, stopped MPI2 and resumed MPI3) plus
  two expected invalid-input refusals. Intermediate and final gas snapshots
  match byte-for-byte across serial/MPI and changed-rank restart. Evidence:
  `application-runs/reacting_restart_wp1ly2_p/checks.json`.
- Python syntax, `git diff --check` and `tools/check_project_contract.py --source .`
  passed. `verification.json` records commands and final source/binary hashes.

The fixed-reference comparison completed all three root steps, but **fails
20 unchanged field limits**. Evidence:
`reference-runs/matched_eb_amr_xpqz290o/checks.json`. At 300 ns the maximum
absolute differences are density `1.238252266e-7 kg/m3`, x-momentum
`6.221881769e-4 kg/m2/s`, total energy `0.787347589 J/m3`, held energy
`0.630854270 J/m3`, and temperature `2.572305657e-3 K`. This is a modest
improvement over the preceding comparison, not trajectory parity. The maximum
mass/element inventory difference from the reference is `5.99e-14` relative;
the energy-budget residual is at most `1.31e-16`, and the H2-consumption
relative difference is `2.24e-5`. No tolerances or reference hashes changed.

## Remaining diffusion discrepancy

A separate initial-rate diagnostic subtracts the pinned native transport-off
MOL rate from its full-transport rate. Valid-cell input density, momentum and
total energy are identical between those native evaluations. The prepared
Fortran transport rate still differs, particularly at fine-patch boundaries.
`transport/comparison.json` records the maximum errors and their cells.

The pinned `Source/Diffusion.cpp` interpolates the sum of hydro and diffusion
fluxes from face centers to open-face centroids. The Fortran molecular
transport kernel currently omits that interpolation. A copied kernel with
centroid interpolation reduces the coarse-level maximum x-momentum rate
error from `325.53` to `5.59e-4 kg/m2/s2`. It also reduces the fine-level
energy-rate error from `3.19e5` to `3.50e3 J/m3/s`, while the fine z-momentum
rate error remains about `2.99e3 kg/m2/s2`, largest at its upper interface.
Evidence: `centroid/comparison.json`, copied source, commands and hashes.

That diagnostic is **not production code**: it interpolates after the existing
species limiter and therefore does not establish positivity of the corrected
fluxes. The next implementation must preserve the species/enthalpy limiting
contract while adding centroid transport fluxes, then isolate the remaining
fine-interface viscous stencil/halo difference. Primitive energy substitution
alone did not improve the coupled trajectory and is not the selected fix.
The full configuration matrix, installed-binary audit, convergence and physical
validation were not run in this increment.
