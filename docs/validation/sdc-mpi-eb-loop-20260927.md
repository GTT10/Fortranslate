# Common-step MPI EB tree SDC — 2026-09-27

## Source and behavior

Implementation `d26217c1`, based on `39cb0863` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`advance_mpi_eb_tree_sdc_3d` couples the owner-local SDC reaction, molecular
transport and CTU hydro stages. The private initial snapshot is restricted;
intermediate stages preserve provisional parents. Each iteration evaluates
hydro from the initial gas using old/new transport and force rates plus the
reaction estimate, recovers the provisional state, refreshes source rates,
and integrates forced chemistry from the initial state. Final hydro and
trapezoidal transport fluxes are synchronized and refluxed before publication.

Chemical, external-force, exposed-wall and physical-boundary inventories remain
separate. Reaction estimates, state, temperature and all inventories publish
only after collective success. Rank-dependent iteration counts, operator flags,
reaction-estimate presence, timestep, tolerances and external force are rejected.
Covered reaction estimates retain caller values. No small-cell redistribution
or selected-level stepping is added by this common-step driver.

## Completed focused checks

GNU Fortran 13.3.0, configured `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- `test_mpi_eb_tree_full_3d` built successfully.
- CTest `unit_mpi_eb_tree_full_3d_2`: 1/1 passed, 67.51 s test / 67.56 s total.
- New comparison uses three levels, five patches, two ranks, 1000 K initial
  fluid, dt 1e-10 s, two SDC iterations, chemistry, all four transport flags,
  outflow faces and force `[2e4,-3e4,4e4]`. Both round-robin owners and all
  patches on rank 0 (rank 1 idle) were exercised.
- State, temperature, reaction estimates and limiter match serial exactly.
  Separate inventories agree within the asserted relative bound 2e-12.
  Composite conservation residual is `4.2614898476812615e-15` in both layouts
  (asserted bound 3e-14).
- Invalid artificial-viscosity coefficient rejects after initial transport;
  state, temperature, reaction estimates and all five scalar/vector outputs
  retain their input values. A separate rank-dependent force, including on an
  idle rank, rejects before stepping and preserves the same outputs.
- Existing MPI reaction-stage and full split-step tests run in the executable.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_full_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs and captured test output: `build/sdc-mpi-eb-loop-checks/`.
No selected build, MPI EB SDC application/restart/regrid, recursive stepping,
operator-disabled SDC case, full configuration matrix, installed-binary audit,
external solver parity, temporal convergence, scaling or physical validation
was run. This is bounded serial/MPI implementation agreement.
