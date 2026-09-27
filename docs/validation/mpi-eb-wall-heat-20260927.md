# MPI embedded-wall heat donors and composite ledger — 2026-09-27

Implementation: `31e6356f`. `compute_mpi_eb_wall_heat_3d` connects the
[wall heat source](eb-wall-heat-source-20260927.md) to the existing MPI AMR
halo exchange. It requires replicated geometry with two target-resolution
layers, validates its interior against the physical patch, and exchanges
owner-local conserved states and temperatures with matching volume fractions.
The source operates only on the patch interior. Existing composite integration
excludes coarse cells covered by children when summing wall heat.

The MPI SDC transport entry now optionally adds the source and returns integrated
wall energy and a global frozen-coefficient frequency bound. The four arguments
`embedded_wall_temperature`, `extended_geometry`, `embedded_wall_energy` and
`embedded_frequency` must be supplied together. Temperature/enable settings
must agree across ranks and conduction must be enabled. All outputs publish
after collective success. Target depth and ancestor temporal inputs are forwarded
to the existing exchange; recursive wall-heat execution is not yet verified.

Application settings, global time-integrator wall-ledger accumulation, timestep
selection and checkpoint identity remain to be connected. This is not a full
Spray-EB run or reference trajectory comparison.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug configuration
`build/heptane-simple-mpi-debug`, spray ON, tests/SUNDIALS OFF. MPI support was
rebuilt successfully without changing toolchain/backend settings. Build logs:
`build/mpi-eb-wall-heat-build.log` and `build/mpi-eb-wall-heat-final-build.log`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_support --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -J build/mpi-eb-wall-heat-unit -I build/mpi-eb-wall-heat-unit -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules src/parallel/mpi_eb_wall_heat_3d_mod.F90 src/parallel/mpi_eb_tree_transport_3d_mod.F90 tests/unit/test_mpi_eb_wall_heat_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/mpi-eb-wall-heat-unit/check
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/mpi-eb-wall-heat-unit/check
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffree-line-length-none -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules tests/unit/test_mpi_eb_tree_transport_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/mpi-eb-wall-heat-unit/existing-check
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/mpi-eb-wall-heat-unit/existing-check
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Two focused MPI executables passed:

- New test, log `build/mpi-eb-wall-heat-unit/corrected-run.log`: two fine patches
  cover the root. A cut cell in the last x-column of the left patch requires
  temperature donors from the right patch. Distributed and all-on-rank-zero
  layouts (including an empty rank) match a whole-fine-grid source reference:
  `-42.739859567646036 W`, relative gate `1e-12`. Covered coarse wall cells are
  excluded. SDC source routing, integrated wall energy and collective invalid
  temperature rejection pass.
- Existing 3D MPI EB transport test, log `build/mpi-eb-wall-heat-unit/existing.log`:
  default transport, SDC checks and serial/MPI comparisons pass unchanged gates.

The first test compile omitted the MPI module include directory and failed;
adding `mpi_modules` fixed the command. The initial executable then rejected
the wide-halo request because the required `ghost_volume_fraction` output was
absent. The implementation now supplies it. Original failure log `run.log` and
diagnostic `debug.log` are retained; neither is counted as a pass. Linker
executable-stack and MPI X authorization warnings remain; corrected executions
exit zero. No full CTest matrix, installed-binary audit or physical validation
was performed. No Python source changed. The separate full single-jet process
remained live and its binary was untouched.
