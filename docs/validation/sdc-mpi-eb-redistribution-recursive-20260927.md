# MPI recursive SDC redistribution — 2026-09-27

Implementation `a401c241`, based on `75c45f63` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

MPI recursive SDC now accepts `redistribution='flux_redist'`. The selector is
agreed collectively before branching. Owner-local interval registers follow
the serial schedule: a transfer from level s to level r belongs to min(s,r),
and is applied when that coarser interval closes. Final reflux increments are
redistributed; changes directed to still-open coarser levels remain deferred.
Temperature recovery and fluid-volume restriction follow each closed interval.
All candidates/registers stay private until the root step succeeds.

The MPI application now permits `subcycling=.true.` with flux redistribution.
The existing common-step path and default `none` selector remain available.
The source is shared with selected-mechanism builds; only fixed builds were
checked in this increment.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort`, `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`; existing compiler warnings remain.

- Mechanism provenance, MPI test and fixed serial/MPI application targets built successfully.
- `unit_mpi_eb_tree_sdc_redist_recursive_3d_2`: 1/1 passed,
  43.34 s test / 43.39 s total. Three levels, five patches, reactions, all
  transport operators, external force, two iterations and root dt 1e-10 s.
- Round-robin ownership and all patches on rank zero match serial gas,
  temperature and reaction estimates exactly. Level counts are `[1,2,4]`,
  limiter is exact and ledgers agree within 2e-12 relative comparison.
  MPI composite residual is 4.4427037072052231e-15 for both layouts.
- A rejected hydro stage preserves gas, temperature, reaction estimates,
  ledgers, limiter and step counts.
- Four application runs completed with `sdc_flux_redist.nml` plus subcycling:
  serial continuous, MPI two-rank continuous, MPI two-rank checkpoint after
  one step, and one-rank continuation. Two root steps finish at 2e-10 s.
  Serial/MPI state CSVs are byte-identical. Continuous/restarted MPI state and
  budget CSVs are byte-identical. Maximum relative budget residual across
  checked runs is 9.416036854481e-16.
- Python syntax, project contract and `git diff --check` passed. No failures,
  skipped tests or interrupted runs were counted.

Artifacts: `build/sdc-mpi-eb-redist-recursive-checks/`, including `provenance.log`,
`build.log`, `ctest.log`, `detail.log`, `contract.log`, `run.py`, `app.log`,
`summary.json`, and isolated inputs/logs/CSVs/checkpoint.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_mpi_eb_tree_full_3d pelef_mpi_reactive_eb_tree_3d pelef_reactive_eb_tree_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_mpi_eb_tree_sdc_redist_recursive_3d_2$' --output-on-failure
python3 build/sdc-mpi-eb-redist-recursive-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-mpi-eb-redist-recursive-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

This is implementation, serial/MPI comparison and restart evidence. It does not
establish temporal convergence, increased stable timestep or upstream AMR
redistribution scheduling parity. Intermediate reflux admissibility checks remain.
No selected build, regrid trajectory, large-tree scaling, installed audit, full
matrix or physical validation was run. Tree state redistribution and spray SDC
remain open.
