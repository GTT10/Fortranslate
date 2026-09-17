# Owner-local MPI EB initialization — 2026-09-18

Implementation: `ee77408ee4029dcb17c61450a224121f165aa0d7`, based on
`a4c171c`. Adds `initialize_mpi_eb_tree_3d` to fixed and selected MPI runtimes
and replaces full-tree construction/scatter in the shared MPI application.

## Implemented behavior

The initializer builds and collectively validates replicated patch geometry,
then allocates U/T only for locally owned nodes. The root owner initializes its
patch using the serial ratio-one PCM/EOS kernel. Each child owner requests only
the required parent section through the existing sparse cell exchange, copies
the section's inactive payloads, and uses the same serial PCM kernel for active
cells. Parent-before-child ordering retains repeated EOS recovery at each depth.
Requested cell counts are checked before integer multiplication/allocation.

New geometry and U/T publish together only after collective success. The sparse
exchange requires fetched parent payloads to be finite, including inactive
cells; this is the existing cell-exchange contract. No full refined state tree
or root scatter buffer is allocated. Application tree nodes now hold geometry
metadata only. Initial composite inventories use MPI reduction. Restart topology
activation uses the same owner-local initializer before owned records are read.

The application still constructs the root seed U/T on every rank and keeps it
through restart topology setup, then releases it. Geometry remains replicated.
Thus this removes full-tree state replication, while root seed replication and
the legacy checkpoint boundary scan remain separate limitations.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`, Open MPI wrapper
`/usr/bin/mpifort` and launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-initialize-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-initialize-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-initialize-debug -R '^(regression_mpi_eb_tree_scheduled_3d|regression_selected_mpi_eb_tree_3d)$' --output-on-failure --output-log build/mpi-eb-tree-initialize-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips, total **53.53 s**:

- Fixed explicit scheduled case: **26.45 s**, including 5 -> 1 -> 4 nodes and
  restart on empty ranks at the root-only stage. Maximum normalized inventory
  residual: **7.05347627465518679e-15**.
- Selected implicit H2/O2 case: **27.04 s**, five-node tree. Maximum normalized
  inventory residual: **3.32177427830221154e-15**.

Serial/two/three-rank final state CSVs and checkpoint patch records match exactly.
Same-rank restart budgets match exactly; changed-rank and serial/MPI restart
directions pass unchanged 3e-14 inventory-scaled budget bounds and independent
fluid-volume integrals. Existing corrupt state, terminator, policy and selected
identity/composition rejection checks preserve sentinel outputs.

Provenance, configure/build, Python syntax, project contract and whitespace
checks pass. No failed build/test attempt occurred; compiler warnings remain.
No Python source changed. These are focused application checks, not separate
initializer fault-injection, peak-memory measurement, install, scaling, full-matrix
or physical-validation evidence.

Logs: `build/mpi-eb-tree-initialize-checks/{configure,build,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-initialize-debug/tests/{mpi_eb_tree_scheduled_3d,selected_mpi_eb_tree_3d}/`.
