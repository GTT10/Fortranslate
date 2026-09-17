# Owner-only MPI EB root seed — 2026-09-18

Implementation: `0465e5e390cbbb5d13219961370fc22f0f810c98`, based on
`1cae205`. The shared fixed/selected MPI EB tree application now partitions
before constructing its analytic root seed. Only the root patch owner allocates
and initializes that seed; other ranks pass zero-sized arrays to the collective
tree initializer. Seed allocation/initialization success is reduced collectively
before entering the initializer. Local seed arrays are released on return.
Restart topology activation repeats this procedure with the restored layout's
owner map. Serial initialization and the numerical kernels are unchanged.

## Completed verification

Fresh GNU Fortran 13.3.0 Debug/Ninja configuration, `/usr/bin/gfortran`,
Open MPI `/usr/bin/mpifort` and `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-seed-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-seed-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-seed-debug -R '^(regression_mpi_eb_tree_scheduled_3d|regression_selected_mpi_eb_tree_3d)$' --output-on-failure --output-log build/mpi-eb-tree-seed-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips, **55.53 s** total:

- Fixed explicit scheduled tree: **25.36 s**, including 5 -> 1 -> 4 nodes.
- Selected implicit H2/O2 tree: **30.14 s**.

The existing checker compares serial/two/three-rank final state CSVs and
checkpoint patch records exactly, same-rank restart budgets exactly, and
changed-rank budgets within the unchanged inventory-scaled tolerance.
It also checks serial/MPI restart in both directions and malformed checkpoint
rejection with protected output sentinels.

Provenance, configure/build, Python syntax, project contract and whitespace
checks passed. A repeat build after the API comment edit reported no work to do.
No failed build/test attempt occurred; compiler warnings remain. No separate
allocation-failure injection, peak-memory measurement, scaling, install or
physical-validation claim is made. Geometry remains replicated and rank zero
still scans legacy checkpoint record boundaries.

Logs: `build/mpi-eb-tree-seed-checks/{configure,build,rebuild,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-seed-debug/tests/{mpi_eb_tree_scheduled_3d,selected_mpi_eb_tree_3d}/`.
