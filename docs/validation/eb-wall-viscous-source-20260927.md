# Stationary embedded-wall momentum source — 2026-09-27

Implementation: `ac2bc929`. `compute_eb_wall_transport_source_3d` now combines
independently enabled conduction and stationary no-slip traction using the same
geometry, recovered primitive variables and active donor stencil. The former
heat-only API remains as a wrapper. Simple transport supplies its existing bulk
viscosity; the analytic transport route retains zero bulk viscosity. Normal
traction uses `4*mu/3+bulk`, tangential traction uses `mu`. Source rates divide
extensive force by fluid volume. No mass/species transfer or mechanical wall
energy is introduced for the stationary wall. The frozen frequency bound also
includes the viscous row sum divided by fluid-volume mass.

`compute_mpi_eb_wall_transport_3d` exchanges both velocity and temperature donors
and returns the composite force/heat ledger. Its heat-only wrapper preserves
the existing SDC/app path. Covered coarse cells are excluded through the existing
composite integrator. Full traction time integration and application no-slip
controls are still pending; the new force-source API alone does not establish
evolved viscous-wall parity.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug build
`build/heptane-simple-mpi-debug`, MPI/spray ON, tests/SUNDIALS OFF.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_support --parallel 2 > build/eb-viscous-source-build.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules tests/unit/test_mpi_eb_wall_heat_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/eb-viscous-source-unit/check
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/eb-viscous-source-unit/check > build/eb-viscous-source-unit/run.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Build and commands pass. The focused MPI executable adds a 2 m/s tangential
flow to the prior two-patch plane fixture. The force is
`[0,-5.9822404766038384e-4,0] N` for both distributed owners and all-on-rank-zero
layouts. It matches a whole-fine-grid evaluation within `1e-12` relative,
excludes covered coarse contributions, and returns zero mass/species/energy
exchange. This fixture uses analytic transport; fitted bulk traction retains
the separately checked kernel and coefficient implementation.

The same executable retains cross-owner heat, empty-rank, collective failure
and evolved common/recursive heat checks; their cooling and inventory residual
remain `-4.2739859565118277e-6 J` and `1.7922725687301810e-15`.
The standalone heat-source unit also passed after the shared implementation
change (`build/eb-wall-viscous-unit/heat-check`), including analytic heat, SDC
routing and late-donor rollback. This was compiled directly from the changed
kernel/source modules and the preceding core archive before the full core rebuild.

No tolerance changed, Python source changed, full CTest matrix, selected-runtime
application rebuild, full reference-solver comparison or physical validation
was performed in this increment. The existing long decane jet remained live.
