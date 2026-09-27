# Rotating-frame source and MOL flux kernels — 2026-09-27

## Implemented boundary

Source `f38f1a75ad2abecf4d51f29392d932b7981696c6`, based on `43bf42cc`, adds
constant, axis-aligned rotation calculations and optional rotating-energy
correction to the existing MOL interface flux. **Application/time-integration
support is not yet connected.** The pinned reference is PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; `Setup.cpp` requires 3D and MOL
when rotation is enabled.

`rotating_frame_3d_mod` provides a frame (Fortran axis 1/2/3, angular speed in
rad/s and a point on the axis), centrifugal potential at a supplied cell-center
position, relative-momentum source, and reversible energy conversion:

```
Phi = |Omega x (x-axis_point)|^2 / 2
S_m = -rho [2 Omega x v_rel + Omega x (Omega x (x-axis_point))]
E_rot = rho e + rho |v_rel|^2 / 2 - rho Phi
```

The rotation source has zero mass/species/rotating-energy terms. EOS energy
recovery adds `rho Phi` before subtracting relative kinetic energy. Valid
rotating total energy can be negative. Position/axis inputs are validated and
failed calls retain outputs. Potential evaluation directly sums the two
perpendicular squared distances, avoiding cancellation of a large axial
component; this is algebraically equivalent to the reference radius formula.

`mol_riemann_flux_3d` accepts two optional cell-center potentials. The acoustic
branch selects potential by the Riemann contact-speed sign, selecting the right
side at zero contact speed as in `Riemann.H`. The LF branch includes both
potentials in advective energy flux and the dissipative energy jump. Omitted
potentials retain the existing path. Neither branch alters mass/species flux.

## Focused checks

GNU Fortran and C++ 13.3.0. The existing Debug fixed core was rebuilt in
`build/heptane-simple-mpi-debug` (tests disabled); direct test programs use
`-g -fcheck=all -ffree-line-length-none`.

- `test_rotating_frame_3d`: 61 checks pass, covering all axes and rotation
  signs, centrifugal/Coriolis directions, zero Coriolis work, shifted origins,
  zero-rotation identity, energy round trips including negative rotating
  energy, and invalid-input rollback.
- `test_mol_riemann_3d`: 82 checks pass, including previous nonrotating cases,
  first/second-order reconstruction in every coordinate direction, analytic
  rotating acoustic/LF energy flux and invalid-potential rollback.
- Existing `test_reactive_pelec_riemann` passes after the primitive-flux API
  exposes its already-computed contact speed.
- `check_pinned_rotating_mol.py` compiles unchanged `riemann` and
  `laxfriedrich_flux` bodies obtained with `git show` at the pinned revision.
  Sixty cases/720 flux values pass the existing Riemann scaled-error gate
  `3e-12`; worst error is `4.362242981416919e-15`. Cases cover both fluxes,
  positive/negative/zero velocities, supersonic flow, unequal thermodynamic
  states and zero/unequal centrifugal potentials. The EOS adapter uses exported
  seven-species NASA7 coefficients in SI, as in the existing Riemann oracle.
  This compares interface arithmetic, not the full upstream EOS or a trajectory.
- Selected serial EB application rebuild succeeds (build only). Provenance,
  project contract, touched Python syntax and `git diff --check` pass.

No CTest matrix, installed-binary audit or rotating CFD/application test was
run. Existing compiler/executable-stack warnings remain.

## Reference and artifacts

Pinned source SHA-256:

| Source | SHA-256 |
| --- | --- |
| `Source/RotSource.cpp` | `4def1b018a8b3ead4dd794ec4c271386d672ca70436bf7d7018b1ab77c0b6cd8` |
| `Source/PeleCUtilities.H` | `2a77e64ce3cf3cb4eb00716ae6de66b9e470f007a040c0c812878142ea6e940d` |
| `Source/Riemann.H` | `073ae14b2f7925de2858ef312d02fc50f033c2a344dccc5b24ecde1f77ed8253` |
| `Source/Setup.cpp` | `2d072a2bd261e871bf57556a1ad41a1e41612627bced999bc8d143678f6c760d` |

`build/rotating-frame/commands.txt` records commands. `check.log`, `mol.log`,
`default.log`, `build.log` and `selected-build.log` retain outcomes.
`reference/` contains the source extraction, input states, exported thermo,
compiled adapter, reference/Fortran outputs and `report.json`.

| Executable | SHA-256 |
| --- | --- |
| Rotation unit | `31aa0db68987de073cd5c86646f6151c6495b5b35fc82374bad1804c1357f7c2` |
| MOL unit/stream | `55af6a6ef8197ff7fe15fc83490f14f5bdea82893cf28b49a0091e39dce1237a` |
| Existing Riemann unit | `f3312324e1f61ff409cf1a9a85a0d7f134ea458b4992c0da8183145378ce2a2b` |
| Selected serial app, build only | `23da2a533c1e21e759392bc372a53fae98bab9177f802ef4005f15f5853dcc2b` |
