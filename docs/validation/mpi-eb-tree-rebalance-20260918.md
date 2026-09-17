# Measured-work MPI EB tree rebalancing — 2026-09-18

Implementation: `3092113996f3007757aa79972d1d62810194f609`, based on
`7ff1cd0`. Fixed/selected MPI EB applications expose `rebalance_interval_steps`
(default 0) and `rebalance_minimum_gain` (default 0.1). Owner-local reaction,
hydro and transport kernels use the existing process-local work sampler,
including recursive hydro/transport substeps. Each retry starts a fresh sample;
only an accepted step can trigger rebalancing. Allocation success is collective.

The existing largest-work-first rebalancer evaluates the measured per-patch
costs and migrates U/T together only for positive estimated peak-load improvement
meeting the gain threshold. Geometry remains replicated. Sampling stops before
rebalancing and regridding; a later sample uses the new layout. The controls are
performance policy rather than checkpoint identity, so restarted runs can
choose different owners. A nonzero interval is rejected by the serial executable.

Timing excludes communication, synchronization, wall-source bookkeeping and
migration. The sampler is process-local, non-reentrant and not thread-safe,
as in the regular MPI tree application. No speedup or scaling claim is made.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja configuration, `/usr/bin/gfortran`,
Open MPI `/usr/bin/mpifort` and `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-rebalance-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-rebalance-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-rebalance-debug -R '^(regression_mpi_eb_tree_rebalance_3d|regression_selected_mpi_eb_tree_rebalance_3d)$' --output-on-failure --output-log build/mpi-eb-tree-rebalance-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips, **107.20 s** total:

- Fixed explicit, recursive subcycling, scheduled 5 -> 1 -> 4 nodes:
  **69.30 s**. Both uninterrupted MPI2 and MPI3 logs record two migrations;
  same/changed-rank and scheduled-event restarts also record migrations.
- Selected implicit H2/O2, common stepping: **37.86 s**. MPI3, changed-rank
  restart and serial-to-MPI restart each record one migration.

The tests request an interval of one and zero minimum gain. Timing-dependent
migration counts above are observations from this run, not deterministic test
requirements. Final state CSVs and checkpoint patch payloads match serial bytes
exactly. With rebalancing enabled, same-rank restart budgets use the existing
3e-14 inventory-scaled comparison instead of byte equality because fresh timings
can change MPI reduction grouping. State byte equality, independent fluid-volume
integrals, the 5e-10 conservation bound and all malformed-checkpoint checks remain
unchanged. The original no-rebalance tests retain budget byte equality.

Provenance, configure/build, Python syntax, project contract and whitespace
checks passed. After test-registration edits, a second build regenerated CMake
and reported no compilation work. No build/test attempt failed; compiler warnings
remain. Retry sampling lifetime was inspected in code but no separate retry or
allocation failure was injected. This is focused application evidence, not a
full matrix, install audit, timing benchmark or physical validation.

Logs: `build/mpi-eb-tree-rebalance-checks/{configure,build,rebuild,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-rebalance-debug/tests/{mpi_eb_tree_rebalance_3d,selected_mpi_eb_tree_rebalance_3d}/`.
