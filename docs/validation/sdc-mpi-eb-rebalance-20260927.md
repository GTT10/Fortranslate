# MPI EB SDC source rebalancing — 2026-09-27

Implementation `c1096a66`, based on `e44d9dcd` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`rebalance_mpi_amr_tree_with_field_3d` packs an additional signed cell field
with conserved gas values and reuses the existing whole-patch transfer. Gas,
temperature, field and owners publish together after collective success. There
is no interpolation or gather. The MPI EB SDC application supplies its reaction
estimates and updates the checkpoint record after successful rebalancing.
Chemistry-disabled SDC continues through the gas-only rebalance route.

The existing rebalance selector also now exits when its candidate owner map
is unchanged. Previously differing summation order could report a tiny positive
gain and a rebalance event even on a one-rank communicator.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance, required builds, Python syntax, project contract and
  `git diff --check` passed.
- Final CTest `unit_mpi_amr_patch_tree_distribution_3d_2`: 1/1 passed,
  0.49 s test / 0.54 s total. The earlier run also passed.
- A three-patch, three-level synthetic fixture starts all work on rank 0 and
  uses costs `[9,6,3]`. Gas, temperature and negative extra-field values move
  exactly to the expected owners `[0,1,1]`. Non-owned fields remain unallocated.
- Gain thresholds, invalid owner costs, rank-dependent thresholds and invalid
  temperature retain the existing rejection behavior. Failed migration retains
  signed source values; unchanged/zero work retains owners.
- Three initial application invocations enabled sampling every step with zero
  minimum gain: two-rank continuous, two-rank checkpoint stop, one-rank resume.
  State and budget CSVs were byte-identical across continuous/resumed execution;
  the state also matched the prior serial app fixture. No actual migration was
  observed in the two-rank application runs.
- The initial one-rank resume reported an unchanged-owner event. After fixing
  that selector condition, the unit test and one-rank resume were rerun against
  the final build. The event disappeared; state and budget still matched the
  previously completed continuous run exactly. Other application invocations
  were not repeated after this guard was added.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_distribution_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_distribution_3d_2$' --output-on-failure
python3 build/sdc-mpi-eb-rebalance-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-mpi-eb-rebalance-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The final isolated resume used the same runner's `run` function with
`resume1_final`, one rank and `stop2/saved.chk`, then compared both CSVs and
asserted absence of a rebalance event. Logs, generated inputs, checkpoint,
runner and CSVs: `build/sdc-mpi-eb-rebalance-checks/`. `summary.json` describes
the initial runs, including the subsequently diagnosed false event;
`final-resume.log` and `ctest-final.log` record the final checks.

Actual field migration is evidenced by the unit fixture, not by these measured
application timings. No selected build, full configuration matrix, installed
binary audit, SDC regrid, recursive stepping or physical validation was run.
