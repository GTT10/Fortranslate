# MPI composite EB flux redistribution — 2026-09-27

Implementation `ad5735ef`, based on `4bc6ca55` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`redistribute_sdc_mpi_eb_tree_flux_rate_3d` implements the serial composite
flux-rate operator with owner-local fields. Geometry and leaf masks remain
replicated. Two sparse cell exchanges fetch neighbor rates for owned cut cells,
then fetch donor excess for owned recipients. No full gas or rate tree is
gathered. Excess accumulation preserves serial donor order.

Peer, coarse/fine and paired periodic neighbors use physical fluid-volume
weights. Hidden parent cells do not exchange excess. Solid rates are zeroed;
unused solid NaNs are ignored. Metadata, geometry, boundary controls, field
shapes and fluid values are validated collectively. Outputs publish only after
all ranks succeed, including ranks with no owned patches.

Common-step MPI SDC now accepts `redistribution='flux_redist'`, collectively
agrees the selector, and redistributes old/new transport rates, hydro rates
and the final reflux increment before temperature recovery and restriction.
The MPI application exposes the existing `sdc_flux_redist.nml` fixture. Its
checkpoint records preserve redistribution identity through the existing
embedded patch serializer. Recursive/selected-level redistribution remains
explicitly unsupported; the default remains `none`.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort`, `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance, MPI support, two test targets and fixed serial/MPI
  application targets built successfully.
- CTest: 2/2 passed, 0.61 s total. Serial operator test 0.02 s;
  `unit_mpi_eb_tree_redistribution_3d_2` 0.54 s.
- MPI operator outputs match serial exactly for root-only, partially refined
  and peer-refined trees, including periodic cases and coarse/fine impulses.
  Both round-robin ownership and all patches on rank zero are exercised.
  The result allocates only owned fields. Rank-dependent boundary labels and
  a nonfinite rate on one owner reject collectively and retain prior output.
- Four application runs completed: serial continuous, MPI two-rank continuous,
  MPI two-rank stop after one step, and one-rank continuation. Three levels,
  five patches, reactions/transport/forcing, two root steps to 2e-10 s.
  Serial/MPI final state CSVs are byte-identical. Continuous/restarted MPI
  state and budget CSVs are byte-identical. Maximum relative budget residual
  across the checked runs is 7.16309624393716e-16.
- Python syntax, project contract and `git diff --check` passed. No unexpected
  build/test failures, skips or interrupted runs were counted.

## Commands and artifacts

Artifacts: `build/sdc-mpi-eb-redist-checks/`, including `provenance.log`,
`build.log`, `build-app.log`, `ctest.log`, `detail.log`, `contract.log`, `run.py`,
`app.log`, `summary.json`, and isolated application inputs/logs/CSVs/checkpoint.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target pelef_mpi_support --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_mpi_eb_tree_redistribution_3d test_eb_tree_redistribution_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target pelef_mpi_reactive_eb_tree_3d pelef_reactive_eb_tree_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_(mpi_)?eb_tree_redistribution_3d(_2)?$' --output-on-failure
python3 build/sdc-mpi-eb-redist-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-mpi-eb-redist-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The selected runtime source list includes the module, but no selected build,
installed audit or full configuration matrix was run. Geometry/neighbor scans
remain replicated and request lists are not deduplicated; large-tree scaling
has not been qualified. Intermediate reflux admissibility checks still apply.
This does not establish a larger stable timestep, temporal convergence, upstream
AMR redistribution parity or physical accuracy. Recursive redistribution,
tree state redistribution and spray SDC remain open.
