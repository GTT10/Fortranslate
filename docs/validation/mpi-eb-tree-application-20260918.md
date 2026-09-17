# Static MPI EB tree application and restart — 2026-09-18

Implementation: `01cb76b7b7d6961bfefd8e7c56dc8e6e75a5d389`, based on
`38c7b81`. `pelef_mpi_reactive_eb_tree_3d` builds from the shared serial source
with `PELEF_MPI_EB_TREE` and links `pelef_mpi_support`; it is included in MPI
install targets. The installed executable itself was not audited here.

## Implemented behavior and boundary

The MPI entry accepts the existing fixed-H2/O2 namelist, geometry, solver,
explicit/implicit chemistry, output policy and checkpoint format for static
trees with fixed common intervals. Existing input and checkpoint validation
remain shared. Automatic intervals, regridding/adaptation and subcycling are
rejected before output. Selected MPI wiring remains open.

Initialization/restart currently constructs or reads a full tree on every rank.
Deterministic whole-patch partitioning then scatters U/T from rank zero. Full
tree state and initialization arrays are released; only owned patch U/T persist
during distributed stepping. Geometry is replicated. Composite balance uses
distributed exposed-cell reductions. Time-step retries preserve all distributed
state through the transactional coupled API.

Checkpoint and final CSV writes gather U/T to rank zero and use the serial
writer. Saved full-tree state is released after checkpoint publication. This
is a root-gathered I/O path, not distributed I/O. The unchanged record format
permits changed-rank restart and serial/MPI interchange. MPI reduction grouping
can change budget roundoff even when state CSVs match exactly. The shared staged
writer retains its existing non-atomic final-copy limitation.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build; compiler `/usr/bin/gfortran`, OpenMPI
wrapper `/usr/bin/mpifort`, launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-app-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-app-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-app-debug -R '^(regression_mpi_eb_tree_application_3d|regression_eb_tree_application_3d)$' --output-on-failure --output-log build/mpi-eb-tree-app-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 tests passed**, zero test failures/skips, total **16.40 s** (serial
regression 2.52 s, MPI regression 13.86 s). The fixed explicit-chemistry fixture
has five nodes, three levels and 272 exposed cells.

- Uninterrupted serial, two-rank and three-rank state CSVs match byte-for-byte.
- Two-rank stop/restart matches both full state and budget CSVs byte-for-byte.
- Two-to-three-rank restart, MPI-to-serial restart and serial-to-three-rank
  restart match the uninterrupted serial state CSV exactly.
- Independent volume integrals agree with final budget inventories within
  3e-14. All budget columns agree across these runs within 3e-14 scaled by
  initial inventory. Maximum normalized residual is
  **3.18097960520324641e-15**, below the existing 5e-10 budget gate.
- Corrupt late checkpoint marker, different interval policy and unsupported
  automatic mode fail without changing sentinel state/budget/checkpoint files.

The initial compilation found a local `mpi_status` name colliding with the
imported `MPI_Status` type. Renaming it `mpi_ierr` fixed the build. No test
failed. Provenance, Python syntax, project contract and whitespace checks pass.
Compiler warnings remain. This is not full-matrix, installed-binary, selected
MPI, large-tree, implicit-application or physical qualification.

Logs: `build/mpi-eb-tree-app-checks/{configure,initial-build,build,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-app-debug/tests/mpi_eb_tree_application_3d/`.
