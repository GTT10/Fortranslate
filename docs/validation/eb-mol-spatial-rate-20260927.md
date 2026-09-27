# Embedded-boundary MOL spatial rate — 2026-09-27

`compute_mol_eb_hydro_3d` now assembles the MOL hydro operator on a prepared
patch with three same-resolution halo cells. It recovers active fluid
primitives, evaluates first- or second-order MOL face-center Riemann fluxes,
and uses the existing reference-ordered EB centroid interpolation and
aperture-weighted divergence. Existing pressure-only embedded-wall sources
complete the conservative rate. Returned face fluxes and boundary/wall ledgers
are dt-integrated; the returned cell rate remains per second.

The operator checks geometry, dimensions, active states and final finiteness
before publishing any output. Covered states are not recovered. Connectivity
for characteristic slopes follows the corresponding face apertures. There is
no CTU time tracing, transverse correction or artificial-viscosity pass.
It uses the existing stationary pressure-wall model; prescribed EB problem
states and rotating frames are not implemented by this routine.

## Completed checks

GNU Fortran 13.3.0, existing Debug core build:

```sh
gfortran -c -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -I build/heptane-simple-mpi-debug/modules -J build/eb-mol-unit src/eb/eb_mol_hydro_3d_mod.F90 -o build/eb-mol-unit/mol.o
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -I build/eb-mol-unit -I build/heptane-simple-mpi-debug/modules -J build/eb-mol-unit tests/unit/test_eb_mol_hydro_3d.F90 build/eb-mol-unit/mol.o build/heptane-simple-mpi-debug/libpelef_core.a -o build/eb-mol-unit/check
build/eb-mol-unit/check
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_core --parallel 2 > build/eb-mol-build.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All passed. One focused executable covers a 4³ central region of an inclined
plane fixture with three halo cells. Both reconstruction orders retain
stationary pressure balance and produce nontrivial rates for a varying
density/pressure/velocity state. The integrated rate equals boundary plus
wall exchange with maximum absolute component residual
`3.4694469519536142e-18`; second-order varying-flow residual is
`1.3234889800848443e-23`. Covered cells contain NaN states to check that they
are not recovered. An invalid active state leaves rates, both ledgers and all
three previously allocated face-flux arrays unchanged.

The unit is registered with CTest and the source is listed in core and selected
EB builds. Only the core target and direct unit executable were run; no full
matrix or independent C++ operator comparison was performed. No Python source
changed. Physical boundary preparation, small-cell redistribution, MPI exchange,
MOL old/new time updates and particle/application coupling remain to connect;
this is not yet a usable MOL application.

The independent jet continuation remained live at about 12 minutes elapsed,
99.9% CPU, with a step-2200 checkpoint. Full completion remains unverified.
