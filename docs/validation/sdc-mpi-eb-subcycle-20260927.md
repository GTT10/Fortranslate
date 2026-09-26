# Recursive MPI EB SDC — 2026-09-27

Implementation `2d5455e8`, based on `5d835e2c` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`advance_mpi_eb_tree_sdc_subcycle_3d` now advances owner-local levels recursively.
Each child uses its refinement ratio worth of substeps and interpolates ancestor
gas/temperature and predictor sources in time. Integrated child fluxes accumulate
before parent/child synchronization. Gas, temperature, reaction estimates,
ledgers and optional per-level step counts publish only after collective success.

The MPI tree application routes `time_integrator='sdc',subcycling=.true.` to
this driver. Its automatic timestep uses a linear refinement factor for coupled
SDC transport intervals; the existing squared factor remains for split transport.
The coupled-timestep selector is agreed collectively. The shared application
source also serves selected mechanisms; this check built the fixed mechanism only.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance passed. The first driver build failed on a Fortran line
  longer than 132 columns; the signature was wrapped before any tests ran.
  The final four requested targets built successfully.
- `unit_mpi_eb_tree_sdc_subcycle_3d_2`: 1/1 passed, 43.15 s test / 43.20 s total.
  Three levels, five patches, two SDC iterations, all operators and external
  forcing, dt 1e-10 s. Both round-robin owners and all patches on rank zero
  match the serial recursive state, temperature and reaction estimates exactly.
  Level counts are `[1,2,4]`; ledgers agree within the existing 2e-12 relative
  comparison. Composite residual is 3.9314133175908747e-15 for both layouts.
  A rejected hydro stage preserves gas, temperature, reaction estimates,
  ledgers, limiter and counts.
- `unit_mpi_eb_tree_transport_3d_2`: 1/1 passed, 13.12 s test / 13.17 s total.
  Existing common/split comparisons now also exercise coupled recursive CFL
  bounds, including thermal walls and inflow, against the serial estimator.
- Four application runs completed: serial continuous, MPI two-rank continuous,
  MPI two-rank stop after one step, and one-rank checkpoint continuation.
  The fixture is `cases/reactive_eb_tree_3d/sdc.nml` with subcycling and automatic
  timestep enabled, two root steps to 2e-10 s. Serial/MPI final state CSVs are
  byte-identical. Continuous/restarted MPI state and budget CSVs are byte-identical.
  Maximum relative budget residual across these runs is 6.036620341508481e-16.
- Python syntax, project contract and `git diff --check` passed. No skipped or
  interrupted tests were counted.

## Commands and artifacts

Artifacts: `build/sdc-mpi-eb-subcycle-checks/`, including `provenance.log`,
`build.log`, `build-transport.log`, `ctest.log`, `recursive-detail.log`,
`ctest-transport.log`, `run.py`, `app.log`, `summary.json`, `contract.log`, and
isolated application inputs, logs, CSVs and checkpoint.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_mpi_eb_tree_full_3d pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_mpi_eb_tree_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_mpi_eb_tree_sdc_subcycle_3d_2$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_mpi_eb_tree_transport_3d_2$' --output-on-failure
python3 build/sdc-mpi-eb-subcycle-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-mpi-eb-subcycle-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

This is focused implementation evidence. The run does not establish temporal
order, physical accuracy, large-tree scaling, installed-binary behavior or a
full configuration matrix. Tree SDC small-cell redistribution and MPI EB spray
SDC remain open. Regrid/rebalance trajectories with recursive SDC were not run
in this check; earlier common-step records remain separate evidence.
