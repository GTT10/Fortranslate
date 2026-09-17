# Indexed MPI EB checkpoint restart — 2026-09-18

Implementation: `bc5947425a6534b52108b76224f3ebf99d7bbac3`, based on
`63b305e`. MPI saves use `PELEF_EB_TREE_CHECKPOINT_V2`, retaining the existing
topology/inventory header and serial patch records. One positive byte length per
node precedes the patch payloads. Rank zero reads these lengths, checks bounded
offset arithmetic, seeks to the final marker and checks EOF; it no longer scans
all patch text. Owners read their own records and verify actual end positions
against the index before publishing state. Existing geometry, inventory,
mechanism identity and averaging checks remain in force.

Serial readers accept V1 and V2 and verify consumed V2 patch lengths. Serial
writers retain V1. MPI readers retain the bounded-memory boundary scan for V1,
so existing serial saves and older MPI files remain supported. V2 is a byte
index: changing patch widths or newline encoding requires rebuilding the index.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`, Open MPI
`/usr/bin/mpifort` and `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-index-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-index-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d test_selected_full_reactive_eb_tree_3d test_selected_full_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-index-debug -R '^(regression_mpi_eb_tree_scheduled_3d|regression_selected_mpi_eb_tree_3d)$' --output-on-failure --output-log build/mpi-eb-tree-index-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-index-debug -R '^regression_mpi_eb_tree_scheduled_3d$' --output-on-failure --output-log build/mpi-eb-tree-index-checks/ctest-fixed-rerun.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

First run: **2/2 passed**, no failures/skips, **67.57 s** total (fixed scheduled
**30.70 s**, selected implicit H2/O2 **36.85 s**). An additional corruption case
was added while the fixed checker was already running, before the selected
checker started. The fixed case was therefore rerun against the final checker:
**1/1 passed**, **32.43 s**, **32.46 s** total. These are two unique regression
tests, not three independent configurations.

Final coverage includes exact serial/two/three-rank final state CSVs and patch
payload bytes, same-rank restart budget equality, unchanged changed-rank budget
bounds, MPI V2 -> serial and serial V1 -> MPI restart, scheduled 5 -> 1 -> 4 node
restart and selected mechanism identity checks. Both serial and MPI reject
negative, oversized and shifted lengths, including paired opposite shifts that
preserve total bytes and the footer position. MPI failure cases preserve result,
budget and checkpoint sentinels. The negative-density patch fixture preserves
its byte length so the index does not mask the payload validation.

Provenance, configure/build, Python syntax, project contract and whitespace
checks passed. The repeat build reported no work to do. No build/test attempt
failed; compiler warnings remain. This is focused restart evidence, not a
scaling benchmark, full configuration matrix, installed-binary audit or physical
validation. Geometry remains replicated; V1 still requires the root scan.

Logs: `build/mpi-eb-tree-index-checks/{configure,build,rebuild,ctest,ctest-fixed-rerun}.log`.
Runtime artifacts: `build/mpi-eb-tree-index-debug/tests/{mpi_eb_tree_scheduled_3d,selected_mpi_eb_tree_3d}/`.
