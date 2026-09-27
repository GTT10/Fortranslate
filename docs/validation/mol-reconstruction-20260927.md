# MOL reconstruction kernel — 2026-09-27

The pinned Spray-EB input uses `do_mol=1`. The source in
`Source/Advance.cpp::do_mol_advance` evaluates old/new spatial rates and combines
them with trapezoidal weights; it does not use CTU characteristic time tracing.
Fortran previously had no MOL-specific reconstruction. The new
`src/hydro/mol_reconstruction_3d_mod.F90` implements the acoustic, transverse
and species-contact slopes from pinned `Source/MOL.H::mol_slope` and the
primitive face construction in `Source/MOL.cpp::pc_compute_hyp_mol_flux`.
The reference remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Its API accepts `rho,u,v,w,p,Y(:)` primitives, central sound speed and explicit
left/right EB connectivity. Inactive donors are not accessed. Slopes use the
reference monotonized centered limiter. The sign test avoids multiplication
overflow by comparing signs directly; this is mathematically equivalent for
ordinary finite differences, not a claim of bitwise equivalence for underflow.
The face routine reconstructs partial species densities, sums density and
normalizes mass fractions. All outputs are transactional, and inadmissible
faces are rejected. The kernel supports gas species; extra passive scalar and
rotating-frame extensions are not included here.

## Completed checks

GNU Fortran 13.3.0, existing Debug core build. Commands:

```sh
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -I build/heptane-simple-mpi-debug/modules -J build/mol-reconstruction-unit src/hydro/mol_reconstruction_3d_mod.F90 tests/unit/test_mol_reconstruction_3d.F90 -o build/mol-reconstruction-unit/check
build/mol-reconstruction-unit/check
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_core --parallel 2 > build/mol-reconstruction-build.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All passed. One focused executable checks 30 assertions: analytic linear
acoustic/contact slopes and both face states in three directions, normalized
species, inactive NaN donors, invalid-active-donor rollback, extrema limiting
and nonpositive-face rollback. The test is also registered with CTest; the
CTest suite was not run. Core and selected-runtime source lists include the
module, but only the core target was built in this increment.

This is a reconstruction kernel, not an implemented MOL application path or
compiled-reference parity result. Riemann flux assembly, EB face/wall fluxes,
old/new source integration, particle staging and application/restart routing
remain to connect. Existing CTU behavior is unchanged. No full matrix or
Python implementation changed. The independent jet continuation remained live
at approximately 4 minutes elapsed; no new checkpoint or completion was yet
observed.
