# Prepared EB molecular transport face fluxes, 2026-10-20

Prepared EB molecular transport now reconstructs fluxes from face centers to
open-face centroids before applying its conservative species limiter. This
closes the missing operation identified in the preceding
[halo connection](eb-mol-transport-halos-20261020.md). It also uses both cells'
tangential velocity stencils at nonperiodic prepared patch interfaces, where
the old boundary path used only the adjacent interior cell.

## Implementation

`src/eb/eb_reactive_transport_3d_mod.F90` uses the existing pinned MOL centroid
stencil for the conservative flux and a separately carried species-enthalpy
flux. Applying the same interpolation weights to both lets the existing
common species-face factor reduce species transport and its enthalpy together;
heat conduction and viscous work are not multiplied by that factor. The
limiter measures outgoing stock from the interpolated flux, before publishing
rates or integrated faces. Impermeable and disabled physical faces retain
their restrictions after reconstruction.

The prepared wrapper computes an extra raw-flux margin and exports the same
corrected grow-3 RHS support and valid integrated faces as before. Its public
shapes, held-energy EOS, and transactional output behavior are unchanged.
The ordinary one-layer transport API retains its face-centered behavior;
centroid reconstruction requires the prepared physical-box neighborhood.
The shared kernel serves serial and MPI callers.

At prepared interfaces, normal gradients retain the correct lower-to-upper
orientation on both boundary sides, and transverse derivatives now average
both cells. Physical wall conditions are still applied after that gradient
construction.

## Independent checks and diagnosis

The existing transport-halo fixture now also checks analytic Fourier flux for
`T = 800 + 100*x*z` at a planar cut-face centroid, and normal/shear stresses
for `v = x*y` at both lower and upper patch interfaces. A large-step binary
species-diffusion case activates the limiter and checks nonnegative species,
closed-box species inventories, zero summed diffusive mass flux and the
species-enthalpy flux relationship. Existing physical-wall, periodic,
held-EOS and late-failure rollback checks remain in the same executable.

Two development attempts failed in the new fixture: its first temperature
profile fell below the thermodynamic range in the far halo, and its first
constant-coefficient selection populated bulk viscosity rather than species
diffusion. The fixture was corrected to use an in-range temperature and the
actual rho-D coefficient; production tolerances were not changed. Logs
`kernel.log` and `kernel-final.log` retain those failures. The corrected
`unit_eb_transport_spatial_halo_3d` passes in 4.00 s (`kernel-final-2.log`).

The initial-rate diagnostic reuses the preceding pinned native full-transport
and transport-off rate captures; those captures were not rerun. Fresh Fortran
rates reduce the fine-level maximum z-momentum rate difference from
`2.986e3` to `1.512e-4 kg/m2/s2`, x-momentum from `2.634e3` to
`8.736e-6`, and y-momentum from `2.635e3` to `2.410e-4`. The coarse x-momentum
rate difference falls from `325.53` to `5.591e-4`. Energy-rate differences
remain (`2.539e3 J/m3/s` coarse, `397.03 J/m3/s` fine). These are spatial
operator diagnostics, not coupled-trajectory or physical qualification.
Evidence: `spatial/comparison.json` and its exact commands.

## Reproduction and scope

Artifacts: `build/eb-transport-centroid-20261020/`. Base: dirty `main` at
`0a0e00ebe9138145fbb0e0e745665d0d768df009`; this increment preserves the
preceding halo/re-reflux work. `before/`, `increment.patch`, and source/hash
records separate this change from the other existing modifications.

Build: `build/reacting-hot-diagnosis-20261019/registration`, GNU Fortran 13.3.0
Debug with bounds/runtime and invalid/zero/overflow checks, Open MPI 4.1.6,
CMake 3.28.3 and SUNDIALS 7.2.0. No compiler, backend or build configuration
changed. MPI commands use `/usr/bin/mpiexec`, `OMP_NUM_THREADS=1` and
`HWLOC_COMPONENTS=-gl,-x86`.

## Completed verification

- Mechanism provenance preflight with regeneration passed:
  `build/cantera-venv/bin/python tools/check_mechanism_provenance.py --bundle mechanisms/lidryer.json --regenerate`.
- Build passed for the transport unit, pinned-reference probe, MPI native
  partial recovery and selected serial/MPI applications (`build.log`).
- Two final focused CTests passed, with no skips: the transport-halo test
  (4.00 s) and `unit_mpi_eb_native_mol_partial_recovery_3d_3` (12.55 s).
- Hot CVODE subcycling and changed-rank restart passed (68.77 s): five
  successful runs plus two expected invalid-input refusals. Intermediate
  and final gas snapshots match byte-for-byte across serial/MPI and 2-to-3
  rank restart. Evidence: `application-runs/reacting_restart_pnwmwkek/checks.json`.
- Python syntax, `git diff --check`, and project contract passed. The final
  source and executable identities and exact commands are in `verification.json`.

The fixed-reference comparison completes all three root steps but still
**fails 19 unchanged field limits**. Evidence:
`reference-runs/matched_eb_amr_9xi8p5to/checks.json`.

| Maximum absolute difference at 300 ns | Before | After |
| --- | ---: | ---: |
| Density (kg/m3) | 1.238252266e-7 | 3.702090356e-8 |
| x-momentum (kg/m2/s) | 6.221881769e-4 | 3.201288544e-5 |
| z-momentum (kg/m2/s) | 7.280314833e-4 | 1.160216988e-5 |
| Total energy (J/m3) | 0.787347589 | 0.457020083 |
| Held energy (J/m3) | 0.630854270 | 0.300697281 |
| Temperature (K) | 2.572305657e-3 | 1.966751177e-3 |

The before column is the preceding recorded run, not a rerun in this increment.
Mass/element inventory differences remain below `5.97e-14` relative, the
energy-budget residual is at most `2.98e-17`, and the H2-consumption relative
difference is `2.24e-5`. Source/reference pins and acceptance limits are unchanged.
The full matrix, installed-binary audit, convergence and physical validation
were not run.

## Remaining energy synchronization difference

A separate copied hydro kernel again substitutes total-minus-kinetic energy
only for primitive conversion. With the corrected transport, the temperature
difference is `1.966770e-3 K` and 19 limits still fail, effectively unchanged
(`energy-primitive/runs/matched_eb_amr_dd955ylp/checks.json`). This substitution
was not adopted in production.

A one-root-step stage dump with full transport shows fine-stage gas states
initially identical and maximum energy differences below `4.8e-3 J/m3`
through the second fine predictor. The completed first root step instead has
an energy difference of `0.16722 J/m3`; held energy and temperature differences
are much smaller (`0.01254 J/m3`, `8.07e-5 K`). Stage comparisons and commands
are in `energy-primitive/stages/`.

The same reference export has a fine-level maximum
`abs(rho_E - rho_e - |rho*u|^2/(2*rho)) = 0.16729 J/m3` at 100 ns,
while Fortran is within `1.46e-11 J/m3`. Both start within `4.35e-11`.
At 200/300 ns the native gap remains `0.16161/0.15653 J/m3`.
`closure-energy-gap.json` records this independent exported-field calculation.
The pinned `Source/PeleC.cpp::reflux()` updates coarse and fine states, while
`post_timestep()` subsequently calls `computeTemp()` on its own level only.
Thus fine-level energy synchronization after parent re-reflux is the next
specific comparison target. The stage evidence and code ordering support this
diagnosis; they do not yet prove that changing reset timing resolves all
remaining field failures. Preserve the held-energy high-cancellation fallback
and conservation/rollback contracts when resolving this policy difference.
