# Embedded-wall viscous traction kernel — 2026-09-27

Implementation: `7eed0db4`. `apply_eb_boundary_viscous_flux_3d` extends the existing quadratic wall stencil
with the pinned normal/tangent construction. The normal traction coefficient
is `4*mu/3 + bulk`; the two tangential coefficients are `mu`. Velocities and
wall velocity are transformed into this frame, differentiated with the stencil,
then transformed back. The result is extensive momentum flux (force in SI).
The active-donor threshold remains `1e-14`; invalid active donors or negative
transport coefficients reject without changing the output.

This returns traction only. The pinned routine does not supply a moving-wall
energy flux, and this addition does not invent one. Gas source, MPI exchange,
time-integrated wall impulse, timestep controls and application no-slip input
remain to be connected; the existing embedded heat path is separate.

## Completed checks

GNU Fortran/G++ 13.3.0, focused standalone compilation:

```sh
/usr/bin/gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -J build/eb-wall-viscous-unit -I build/eb-wall-viscous-unit src/core/precision_mod.F90 src/eb/eb_boundary_gradient_3d_mod.F90 tests/unit/test_eb_boundary_gradient_3d.F90 -o build/eb-wall-viscous-unit/check
build/eb-wall-viscous-unit/check
python3 build/eb-wall-viscous-unit/compare_reference.py
python3 -m py_compile build/eb-wall-viscous-unit/compare_reference.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

- One Fortran unit executable passed: 24 signed orientations and quadratic
  scalar/vector fields, analytic normal/tangential traction including bulk
  viscosity, covered-donor rollback and negative-bulk rejection. Absolute
  analytic flux gate remains `2e-14`.
- The comparison driver verifies pinned `Source/EB.cpp` bytes against
  PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, extracts the original
  `pc_apply_eb_boundry_visc_flux_stencil`, and compiles it with a serial AMReX
  interface shim and the preceding constructor harness. Twenty-four cases
  give maximum scaled force error `1.4381418748524587e-14`, below the fixed
  `2e-14` gate. Commands and source hash are in
  `build/eb-wall-viscous-unit/comparison.json`; raw outputs are `reference.txt`
  and `candidate.txt`. This is a compiled kernel comparison, not a full solver.
- Mechanism provenance, Python syntax, project contract and whitespace checks pass.

The first reference-harness compile lacked the `AMREX_D_PICK` macro; adding its
3D definition fixed the harness without altering the extracted function.
No full library/application build, CTest matrix, evolved viscous-wall comparison
or physical validation was run in this kernel increment. The live decane jet
binary was left untouched.
