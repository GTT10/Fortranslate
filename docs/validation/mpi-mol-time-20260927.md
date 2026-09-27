# Nonreacting MPI MOL time update — 2026-09-27

The common MPI tree integrator now accepts optional `mol_mode`, default false.
It evaluates the old MOL hydro/transport rate, forms an Euler predictor,
evaluates the new rate and applies their trapezoidal average. This follows
the nonreacting branch of pinned `do_mol_advance`: one old/new correction,
independent of the SDC iteration count. MOL bypasses CTU hydro evaluation.
Flux ledgers use the same weights; pressure-wall force, embedded heat and
viscous traction use their corresponding old/new stage values.

The existing transport, stage redistribution, common-step flux synchronization,
EOS recovery and transactional publication are reused. The six-layer geometry
cache supplies a three-layer MOL view. Rank controls include the mode. Reaction
and custom predictor-builder combinations are rejected because their MOL stage
semantics are not yet connected. The recursive wrapper and particle/application
controls do not expose this mode yet; their existing CTU behavior is unchanged.

## Completed verification

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug MPI/spray build,
tests/SUNDIALS disabled:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_support --parallel 2 > build/mpi-mol-time-build.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules tests/unit/test_mpi_eb_wall_heat_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/mpi-mol-time-unit/check-final
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/mpi-mol-time-unit/check-final > build/mpi-mol-time-unit/run-final.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All passed. The existing nested two-patch/two-rank wall executable now includes
common MOL heat and viscosity cases plus a nonuniform inviscid MOL update.
MOL wall heat is `-4.2739859565118192e-6 J`, with normalized inventory residual
`1.7922725689661285e-15`. Tangential wall impulse is
`-5.9822404766038380e-11 N s`, zero adiabatic wall energy, and residual
`3.1127680828802201e-17`. The nonuniform inviscid case produces a nonzero
update and closes its ledger to `4.8743201100098756e-15`. Existing CTU common
and recursive heat/viscosity checks remain passing in the same executable.

After that run, inspection found that the optional selected-depth
redistribution-transfer path needed an explicitly zero MOL CTU-transfer field.
That initialization was added and the MPI library rebuilt successfully
(`build/mpi-mol-time-build-final.log`). The tested common-step path is unchanged;
selected-depth/recursive MOL evolution has not been qualified. No full CTest
matrix, application rebuild or external C++ time-update comparison was run.
No Python source changed. The executable-stack linker warning remains nonfatal.

The independent Jet continuation remained live at 21 minutes 39 seconds,
99.9% CPU, and had reached a step-2400 checkpoint. Completion is not established.
