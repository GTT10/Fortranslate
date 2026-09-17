# MPI EB subcycled application and automatic intervals — 2026-09-18

Implementation: `af0c8b8ab84e0f611d8bfb1457c9718421ad2847`, based on
`bc2cf02`. The shared MPI EB application now accepts `subcycling=.true.` and
passes the mode into coupled advancement and automatic timestep selection.

`mpi_eb_tree_timestep_3d` now accepts optional `subcycling`. The mode must agree
across ranks. Hydro patch bounds scale by cumulative refinement ratio; diffusion
bounds scale by its square. Hidden parent updates and global donor maxima remain
included. Comparison precedes multiplication to avoid overflowing inactive huge
bounds. The conservative full-interval diffusion bound is retained for split
transport halves. Absent/false mode preserves common stepping.

The application uses the shared serial checkpoint mode envelope. Static,
scheduled and adaptive topology dispatch now reaches the recursive operator;
restart rejects a changed stepping mode. Root-gathered CSV/checkpoint writes and
temporary replicated initialization/restart remain unchanged.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`, Open MPI wrapper
`/usr/bin/mpifort` and launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-subcycle-app-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-subcycle-app-debug --target test_mpi_eb_tree_transport_3d pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-subcycle-app-debug -R '^(unit_mpi_eb_tree_transport_3d_[23]|regression_mpi_eb_tree_subcycle_(scheduled|automatic)_3d)$' --output-on-failure --output-log build/mpi-eb-tree-subcycle-app-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**4/4 passed**, no failures/skips, total **145.84 s**:

- Two/three-rank units: **5.96/5.85 s**. Across three plane cuts, four/five-node
  trees and distributed/single-owner layouts, **144** timestep comparisons
  (72 common, 72 subcycled) match serial exactly. Each mode exercises hydro-only,
  coupled and diffusion-limited bounds. Rank-disagreeing CFL/subcycling controls
  reject without changing the estimate. Existing transport checks also pass.
- Scheduled subcycling application: **56.55 s**. The tree changes **5 -> 1 -> 4**
  nodes at 1.5e-9/2.5e-9 s and completes four steps to 3e-9 s. Restart before and
  after events agrees with uninterrupted execution. Maximum normalized inventory
  residual: **5.94704866236712754e-15**.
- Automatic subcycling application: **77.44 s**. Hydro CFL 1e-5 selects first
  interval **4.1705245814991907e-10 s**, below the 1e-9 s cap. Five steps reach
  2e-9 s, with minimum accepted interval **3.3179025176661483e-10 s**.
  Maximum normalized inventory residual: **3.18097944570051201e-15**.

Both applications use explicit chemistry and full transport. Serial/two/three-rank
final state CSVs match byte-for-byte, including two-to-three-rank and serial/MPI
restart directions. Same-rank restart budgets match byte-for-byte; changed-rank
budgets satisfy the unchanged 3e-14 inventory-scaled bound. Independent fluid-volume
integrals pass. Corrupt records and changed interval/CFL/subcycling policies
reject while preserving sentinel outputs; stepping-mode diagnostics are checked.

Provenance, configure/build, Python syntax, project contract and whitespace
checks pass. No failed build/test attempt occurred; compiler warnings remain.
Adaptive-plus-subcycled and automatic-plus-scheduled combinations were not
independently qualified here. This is not selected MPI, installed-binary, scaling,
full-matrix, temporal-convergence or physical-validation evidence.

Logs: `build/mpi-eb-tree-subcycle-app-checks/{configure,build,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-subcycle-app-debug/tests/mpi_eb_tree_subcycle_{scheduled,automatic}_3d/`.
