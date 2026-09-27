# Recursive MPI MOL connection — 2026-09-27

Implementation commit `1dc44e97` forwards optional `mol_mode` through all four
recursive level-update branches to the common integrator. Its default remains
false. Existing coarse/fine flux synchronization and redistribution transfer
use the MOL fluxes and rates returned by that integrator.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6; existing Debug MPI/spray build with tests
and SUNDIALS disabled. One manually linked, two-rank focused executable passed:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_support --parallel 2 > build/mpi-mol-recursive-build.log 2>&1
mkdir -p build/mpi-mol-recursive-unit
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules tests/unit/test_mpi_eb_wall_heat_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/mpi-mol-recursive-unit/check
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/mpi-mol-recursive-unit/check > build/mpi-mol-recursive-unit/run.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

The nested planar-wall fixture exercises common/recursive CTU and MOL, then a
nonuniform inviscid MOL step followed by recursive MOL with `state_redist`.
Recursive MOL wall heat is `-4.2739859565118200e-6 J`, with normalized inventory
residual `1.7922725689425341e-15`. Adiabatic tangential impulse is
`-5.9822404766038380e-11 N s`, zero stationary-wall energy transfer and residual
`3.1127680828802201e-17`. The recursive redistributed nonuniform step has
residual `2.6554724194323220e-15`. The executable also retains existing donor,
composite-source and collective-rejection checks.

All commands above passed. Executable-stack and X11 authorization warnings
were nonfatal. No Python source changed. No full configuration matrix,
recursive reacting comparison, application rebuild, particle MOL integration
or C++ temporal parity was performed in this change.

## Separate Jet diagnostic result

The diagnostic replay in `build/decane-jet-failure-trace` reached the terminal
Fortran error. GDB recorded accepted time `0.0004287942691061733 s`, accepted
step `2436`, and next injection index `2144`. The replay used the original
binary and retry-run step-2400 checkpoint; `provenance.json` and `trace.log`
record its identity and output. This localizes the failed interval but does
not identify the rejecting numerical operation or repair the Jet case.
