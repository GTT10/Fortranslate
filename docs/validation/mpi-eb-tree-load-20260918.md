# Owner-local MPI EB checkpoint loading — 2026-09-18

Implementation: `9717e32cf499bd8085f43233d08ed7654855634c`, based on
`3dc3cb5`. Adds `read_mpi_eb_tree_checkpoint_3d` to fixed and selected runtimes
and routes application restart through it after restoring topology metadata.

## Implemented behavior

The application reads the small policy/schedule/adaptation envelope from a
formatted stream, selects the stored topology and distributes initial placeholder
state. The new reader validates topology descriptors and ledger metadata. Rank
zero scans patch start/end positions using a fixed-size line buffer and checks
the tree terminator and end-of-file. Only offset vectors are broadcast.

Each owner seeks directly to its patch records and uses the existing serial EB
checkpoint reader. It verifies record boundaries, immutable geometry, mechanism
and composition context, patch diagnostic metadata and state/element inventories.
Distributed composite inventory and private averaging enforce the existing
5e-10 ledger and 2e-12 parent/child consistency bounds. U/T and checkpoint metadata
publish together only after all ranks succeed. No restored patch payload is
gathered, broadcast or parsed on nonowners.

The existing V1 payload and serial checkpoint interchange are retained. Owners
use independent formatted-stream reads; this is not collective MPI file reading.
Rank zero still scans the legacy text payload to discover offsets. Initial
placeholder construction still temporarily builds the initial tree on all ranks
before distribution; this change removes replication of restored state, not
that separate initialization limitation.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`, Open MPI wrapper
`/usr/bin/mpifort` and launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-load-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-load-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-load-debug -R '^(regression_mpi_eb_tree_scheduled_3d|regression_selected_mpi_eb_tree_3d)$' --output-on-failure --output-log build/mpi-eb-tree-load-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips, total **54.37 s**:

- Fixed explicit scheduled application: **25.45 s**, including restart before
  and after 5 -> 1 -> 4 node changes, with empty ranks at the root-only stage.
  Maximum normalized inventory residual: **7.05347627465518679e-15**.
- Selected implicit H2/O2 application: **28.89 s**, five-node tree and mechanism/
  composition-bound restart. Maximum residual: **3.32177427830221154e-15**.

Serial/two/three-rank final state bytes match, including two-to-three-rank and
both serial/MPI restart directions. Same-rank budget bytes match; changed-rank
budgets and independent fluid-volume integrals pass unchanged bounds. Saved patch
records also retain exact serial byte parity. A new rejection case changes the
last patch's final fluid-cell density to -1 while retaining record markers;
owner-side validation rejects it without changing sentinel output/save files.
Existing corrupt terminator, policy and selected identity/composition rejection
checks also pass.

Provenance, configure/build, Python syntax, project contract and whitespace
checks pass. No failed build/test attempt occurred; compiler warnings remain.
No install, adaptive/subcycled restart combination, scaling, full-matrix,
I/O fault injection or physical-validation evidence is claimed.

Logs: `build/mpi-eb-tree-load-checks/{configure,build,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-load-debug/tests/{mpi_eb_tree_scheduled_3d,selected_mpi_eb_tree_3d}/`.
