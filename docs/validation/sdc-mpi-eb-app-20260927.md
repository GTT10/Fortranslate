# MPI EB tree SDC application and restart — 2026-09-27

## Source and behavior

Implementation `3b105d8f`, based on `9132c5c5` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The shared EB tree application routes MPI `time_integrator='sdc'` to the
common-step SDC driver, initializes reaction estimates only on owned patches,
and uses the existing SDC forcing/reaction checkpoint lifecycle. V4 MPI tree
records retain these fields across changed-rank reads. PeleC CTU and explicit
physical boundary faces remain required.

MPI SDC currently rejects subcycling, scheduled/adaptive regridding and load
rebalancing because those paths still need source-state integration. The MPI
EB spray application continues to reject SDC. Serial recursive/regrid routing
is unchanged. Fixed/selected MPI share the source route; only fixed executables
were built and run in this check.

## Completed focused checks

GNU Fortran 13.3.0, configured `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- Fixed serial and MPI EB tree applications built successfully.
- Four application invocations passed using `cases/reactive_eb_tree_3d/sdc.nml`:
  serial continuous, two-rank continuous, two-rank stop/checkpoint after step 1,
  and one-rank continuation to step 2 (2e-10 s).
- The fixture has three levels, five patches, a z-plane EB, dt 1e-10 s,
  chemistry, molecular transport, external forcing and two SDC iterations.
- Serial/MPI final state CSVs are byte-identical.
- Two-rank continuous and one-rank resumed state and budget CSVs are
  byte-identical. MPI composite error is `2.9751533467500573e-16`; the maximum
  normalized residual over serial/MPI/resumed budget CSVs is
  `6.036624067173817e-16` (asserted bound 5e-10).
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_reactive_eb_tree_3d pelef_reactive_eb_tree_3d --parallel 2
python3 build/sdc-mpi-eb-app-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-mpi-eb-app-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The local runner creates isolated inputs, uses `/usr/bin/mpiexec`, removes
DISPLAY/LD_LIBRARY_PATH/LIBRARY_PATH, and sets `HWLOC_COMPONENTS=-x86`,
`OMPI_MCA_hwloc_base_binding_policy=none` and oversubscription permission.
Inputs, runner, checkpoints, logs, CSVs and `summary.json` are in
`build/sdc-mpi-eb-app-checks/`.

No selected build, automatic-timestep trajectory, disabled-operator SDC case,
new rejection-case run, regrid/rebalance, recursive stepping, full matrix,
installed-binary audit, external solver parity or physical validation was run.
