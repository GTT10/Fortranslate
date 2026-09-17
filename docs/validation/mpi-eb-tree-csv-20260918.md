# Owner-local MPI EB CSV output — 2026-09-18

Implementation: `d7fada13eefe9ca5c46f269e84f5f845bda310fc`, based on
`a492fa4`. Adds `write_mpi_eb_tree_csv_3d` to fixed and selected MPI runtimes
and removes final state gathering from the shared application's CSV path.

## Implemented behavior

Each owner validates and measures its exposed fluid-cell rows before opening
the destination. Covered parents and solid cells are skipped. MPI reduction of
per-node byte counts determines node-ordered offsets; each rank then writes its
own rows with `MPI_File_write_at` using a bounded 32-row buffer. Ordering and
`es26.17e3` formatting match serial. Rank zero writes the header and the separate
small budget CSV. No fluid-state gather or full output text buffer is needed.

Topology/ownership/container shape, replicated geometry, path and column header
agreement, and finite emitted values are checked before truncation. Failure in
this preflight preserves existing output. Write counts and collective completion
are checked. As with the existing regular-tree MPI writer, an I/O failure after
opening/truncation is not a transactional file replacement. Checkpoint saving
still gathers state; initialization/restart still temporarily replicates it.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`, Open MPI wrapper
`/usr/bin/mpifort` and launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-csv-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-csv-debug --target test_mpi_eb_tree_csv_3d pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-csv-debug -R '^(unit_mpi_eb_tree_csv_3d_3|regression_mpi_eb_tree_scheduled_3d|regression_selected_mpi_eb_tree_3d)$' --output-on-failure --output-log build/mpi-eb-tree-csv-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-csv-debug --target test_mpi_eb_tree_csv_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-csv-debug -R '^unit_mpi_eb_tree_csv_3d_3$' --output-on-failure --output-log build/mpi-eb-tree-csv-checks/ctest-unit-rerun.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All three required checks pass across the recorded attempts:

- Fixed scheduled application: **passed, 23.13 s**. Serial/two/three-rank state
  CSVs match byte-for-byte, including restart before/after 5 -> 1 -> 4 node
  changes and serial/MPI interchange. Same-rank budget bytes and independent
  fluid-volume inventory checks pass.
- Selected implicit H2/O2 application: **passed, 26.45 s**. Serial/two/three-rank
  state bytes and restart directions match, with identity/composition rejection
  and unchanged budget checks.
- Three-rank CSV unit: initial attempt failed during fixture construction
  (**0.43 s**) because a plane at z=0.5 aligned exactly with a cell face, which
  the existing geometry builder rejects. Moving only the fixture plane to z=0.61
  retained four solid/four fluid cells. The unit-only rerun **passed, 0.47 s**
  (total 0.50 s). It checks active NaN, path/header/geometry mismatch rejection
  with sentinel preservation, an empty root/empty rank, ignored solid NaNs and
  exactly four emitted fluid rows. The first three-test attempt was 2/3 passed,
  total 50.05 s; no application rerun was needed for the fixture-only correction.

Provenance, configure/build, Python syntax, project contract and whitespace
checks pass. Compiler warnings remain. No Python source changed. These are
focused output/API checks, not distributed checkpoints, install, scaling,
full-matrix or physical-validation evidence.

Logs: `build/mpi-eb-tree-csv-checks/` (configure, build, rebuild, first CTest and
unit rerun). Runtime artifacts: `build/mpi-eb-tree-csv-debug/tests/` under
`mpi_eb_tree_scheduled_3d/`, `selected_mpi_eb_tree_3d/` and the unit CSV.
