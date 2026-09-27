# Stationary embedded-wall traction in MPI SDC — 2026-09-27

Implementation: `bd691530`, based on `00b94bf9`. Common and recursive MPI
embedded-boundary SDC now accept optional `embedded_noslip`. The default is
disabled. Enabling it requires viscosity, uses the existing two-layer wall
donor geometry and adds the full wall source once per transport evaluation.
The wall momentum ledger uses the same old/new trapezoidal stage weights as
the transport update. Conduction remains independently selectable; a stationary
adiabatic wall exchanges momentum without mechanical wall energy. Collective
control agreement and transactional publication remain in place.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug build
`build/heptane-simple-mpi-debug`, MPI/spray ON, tests/SUNDIALS OFF.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_support --parallel 2 > build/eb-viscous-integrator-build.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules tests/unit/test_mpi_eb_wall_heat_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/eb-viscous-integrator-unit/check
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/eb-viscous-integrator-unit/check > build/eb-viscous-integrator-unit/run.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All listed final commands passed. The first build failed on a Fortran line
length limit in the new call; splitting that line fixed it. The focused test
executable passed on two ranks, including both common and recursive evolution
of the two-patch fixture with 2 m/s tangential velocity, viscosity enabled and
conduction disabled. Both return tangential wall impulse
`-5.9822404766038380e-11 N s` over `1e-7 s`, zero wall energy, and normalized
total-inventory residual `3.1127680828802201e-17`. Pressure-wall momentum is
included in the existing ledger; the tangential component isolates drag here.

Existing cross-owner/empty-rank source checks, collective rejection and both
heat-evolution checks also passed in this executable. Heat energy and normalized
inventory residual remain `-4.2739859565118277e-6 J` and
`1.7922725687301810e-15`. The linker executable-stack warning and MPI X11
authorization messages did not prevent successful exit.

No Python source changed, full CTest matrix or reference-solver evolution
comparison was run. This increment does not expose no-slip controls in the
application or particle wrappers; that connection, timestep/restart controls,
and full Spray-EB parity remain open. The independent full decane jet run was
still active and was not modified.
