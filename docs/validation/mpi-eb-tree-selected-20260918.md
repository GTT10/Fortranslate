# Selected-mechanism MPI EB tree — 2026-09-18

Implementation: `a4339f4fef2bcb04e1b479f71f5f6093bf29d499`, based on
`2926130`. Adds bundle argument `MPI_REACTIVE_EB_TREE_3D_NAME` and public target
`pelef_mpi_reactive_eb_tree_3d_selected` when MPI and a selected bundle are enabled.

The selected MPI EB runtime compiles the existing EB MPI operators and sparse
tree exchange against the selected serial runtime. It has an isolated Fortran
module directory and does not link fixed `pelef_core` or `pelef_mpi_support`.
The generated application uses both `PELEF_SELECTED_EB_TREE` and
`PELEF_MPI_EB_TREE`; bundle generation remains through the existing pinned
importer and source checks. The target participates in the existing install
list and MPI runtime coherence audit.

The shared application retains selected composition, temperature-range and
element-conservation checks, and binds mechanism SHA, composition and integrator
to checkpoints. Static/scheduled/adaptive and common/recursive paths use the same
MPI operators. Initial setup/restart temporarily replicates state and CSV/save
gathers on the root; distributed I/O remains separate work.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`, Open MPI wrapper
`/usr/bin/mpifort`, launcher `/usr/bin/mpiexec`, configured H2/O2 full bundle.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/selected-mpi-eb-tree-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON -DPELEF_MECHANISM_BUNDLE=mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-mpi-eb-tree-debug --target pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-mpi-eb-tree-debug -R '^regression_selected_mpi_eb_tree_3d$' --output-on-failure --output-log build/selected-mpi-eb-tree-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 regression passed**, no failures/skips, **23.51 s** (total 23.54 s).
The selected full H2/O2 fixture uses H2/O2/N2 composition, implicit chemistry,
full transport and a static five-node, three-level cut-cell tree. Serial and
two/three-rank state CSVs match exactly. Same-rank restart matches state and
budget bytes; changed-rank and both serial/MPI restart directions match state
bytes and meet the unchanged 3e-14 inventory-scaled budget bound. Independent
fluid-volume integrals pass. Maximum normalized inventory residual:
**3.32177427830221154e-15**.

Selected and element contexts are present in saved checkpoints. Altered bundle
SHA and composition are rejected while preserving sentinel state/budget/save
files, as are the existing corrupt-record and stepping-policy failures.
Both public selected targets compile. The public MPI target's Ninja link command
was checked to include `libpelef_selected_mpi_eb_tree_3d_runtime.a` and exclude
`libpelef_core.a` and `libpelef_mpi_support.a`; the exact command is retained in
`link-audit.log`. Build-time MPI runtime dependency audits pass (ABI 40).

Provenance, configure/build, Python syntax, project contract and whitespace
checks pass. No failed build/test attempt occurred; compiler warnings remain.
Runtime regression uses the fixture targets; public targets were compiled and
link-audited. No install audit, selected subcycling/adaptation qualification,
larger-mechanism run, full matrix, scaling or physical validation is claimed.

Logs: `build/selected-mpi-eb-tree-checks/{configure,build,ctest,link-audit}.log`.
Runtime artifacts: `build/selected-mpi-eb-tree-debug/tests/selected_mpi_eb_tree_3d/`.
