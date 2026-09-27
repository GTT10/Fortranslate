# Time-integrated embedded-wall heat — 2026-09-27

Implementation: `201ec51f`. The MPI common and recursive SDC integrators now accept optional
`embedded_wall_temperature`. The common integrator extracts the two-layer heat
geometry from its existing six-layer hydro cache, evaluates heat at the same
old/provisional stages as molecular transport, and adds the final trapezoidal
wall-energy contribution to the existing wall ledger. Recursive stepping
forwards the option on every level interval; its existing composite ledger
excludes covered coarse cells.

Controls agree across ranks, conduction is required, and gas state and all
ledgers retain their previous values on rejected steps. The option remains a
library input: application namelists, timestep policy, particle-wrapper forwarding
and restart identity still need integration. The separate full single-jet run
continues with its original binary.

## Completed evidence

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug build
`build/heptane-simple-mpi-debug`; spray ON, tests/SUNDIALS OFF.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_support --parallel 2 > build/mpi-eb-wall-integrator-build2.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules tests/unit/test_mpi_eb_wall_heat_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/mpi-eb-wall-heat-unit/integrator-check
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/mpi-eb-wall-heat-unit/integrator-check > build/mpi-eb-wall-heat-unit/integrator.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All commands passed. One focused two-rank executable extends the preceding
[cross-owner donor test](mpi-eb-wall-heat-20260927.md) with actual common and
recursive SDC evolution, two SDC iterations over `1e-7 s`, no chemistry,
nitrogen initially at 500 K and a 400 K embedded wall. Two fine patches cover
the entire coarse root; the wall stencil spans their ownership boundary.

| Time integration | Signed wall heat (J) | Maximum normalized inventory residual |
| --- | ---: | ---: |
| Common interval | -4.2739859565118277e-6 | 1.7922725687301810e-15 |
| Recursive intervals | -4.2739859565118277e-6 | 1.7922725687301810e-15 |

The inventory residual subtracts physical boundary, wall, chemistry and external
source changes. Existing gates were retained; the new evolved checks require
residual below `1e-12` and cooling within `1e-5` relative of the initial frozen
heat rate times the short duration. The executable also reruns source routing,
composite coarse exclusion, empty-rank layout and collective invalid-input checks.

The first build failed because the recursive optional argument lacked its
declaration; `build/mpi-eb-wall-integrator-build.log` retains that diagnostic.
After adding the declaration, the build and execution passed. Existing linker
executable-stack and MPI X authorization warnings remain. No full CTest matrix,
installed-binary audit, C++ evolved comparison or physical validation was run.
No Python source changed.
