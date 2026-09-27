# MOL Lax–Friedrichs interface flux — 2026-09-27

Implementation `bbe3c31e` adds optional `lax_friedrichs=.true.` to
`mol_riemann_flux_3d`. The default remains the acoustic mixture Riemann solver.
Reconstruction, positive-face slope protection and momentum-axis permutations
are shared. The new branch evaluates face EOS sound speeds and conservative
states, then follows pinned `Source/Riemann.H::laxfriedrich_flux`:

- Acoustic impedances weight the contact velocity, including pressure jump.
- The reference `smallu=1e-12` relative velocity mask is retained.
- Diffusion speed is absolute contact velocity plus mean face sound speed.
- Species use the LF flux directly; mass flux is their sum.
- Momentum and total energy use the same dissipation speed.

This differs from choosing the maximum one-sided `abs(u)+c` Rusanov speed.
The source is pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Rotating-frame energy and extra advected/linear scalars are not included in
this mixture-state interface. EB spatial/MPI/application selection still needs
forwarding; adding the interface option alone does not enable a CLI LF mode.

## Completed checks

GNU Fortran 13.3.0; existing Debug fixed-core build. Commands:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_core --parallel 2 > build/mol-lf-build-final.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/modules tests/unit/test_mol_riemann_3d.F90 build/heptane-simple-mpi-debug/libpelef_core.a -o build/mol-lf-unit/check-final
build/mol-lf-unit/check-final > build/mol-lf-unit/run-final3.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final checks passed: one executable, 64 assertions across the existing
acoustic and new LF branches. LF checks include analytic physical flux for
uniform/linear fields, equal-pressure stationary contacts with density jump,
pressure-driven contact speed at zero face velocities, species mass closure
and active-NaN rollback in all three directions. Closed outer donors remain
ignored. Existing CTest registration covers the expanded executable; no full
configuration matrix was run.

The first added stationary-contact test chose a density that put its EOS
temperature below the shared NASA range and correctly failed. The final
fixture uses densities 2 and 1.5 kg/m³ at 200 kPa, within that range; the
original failed `build/mol-lf-unit/run.log` remains. The energy-flux gate stays
`2e-13` scaled. This is analytic kernel evidence, not compiled C++ parity or
full-solver qualification. No Python source changed.
