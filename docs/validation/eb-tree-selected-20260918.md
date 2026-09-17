# Selected-mechanism multilevel EB application — 2026-09-18

Implementation: `63a4fcaf605876b314614d053ac81e2023b83706`, based on
`253b4d6`. Executable: `pelef_reactive_eb_tree_3d_selected`.

The fixed and selected executables share one preprocessed application source.
Selected loading reuses mechanism/transport validation, composition resolution
and initial temperature range checks. Generic source element metadata feeds
checkpoint records and the existing elemental chemical-exchange gate. Bundle
SHA-256, initial composition and integrator are passed into every node record.
Fixed builds reject selected composition inputs. Generated mechanism data is
unchanged.

CMake adds the selected tree frontend option and includes tree/checkpoint
modules in the selected EB runtime, linking its existing regular AMR tree
runtime dependency. The selected executable is included in configured bundle
installation targets; no fresh install audit was run. No compiler, MPI or
optional backend was changed. Shared scheduled/adaptive controls are present,
but the selected checks below exercise a static tree only.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`; selected full H2/O2 fixture, implicit chemistry.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_eb_tree_application_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_tree_3d test_selected_full_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_(selected_eb_tree|eb_tree_application)_3d$' --output-on-failure --output-log build/eb-tree-selected-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_selected_eb_tree_3d$' --output-on-failure --output-log build/eb-tree-selected-checks/selected-ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

Initial **2/2 passed**, zero failures/skips, total 7.10 s (fixed 2.47 s,
selected 4.60 s). After adding explicit selected identity rejection checks,
selected **1/1 passed**, 4.89 s (total 4.92 s). No failed attempt occurred.

The five-node, three-level fixture has 272 exposed cells and runs two 1e-9
steps with implicit chemistry, HLLC and molecular transport. Restart yields
byte-identical final state and budget CSVs; independent CSV volume integrals
agree with budgets. Maximum normalized balance residual is
**3.32177427830221154e-15**. Altered checkpoint bundle SHA-256 and changed
input composition produce the corresponding mismatch diagnostics before
output creation. Existing policy/mode/corruption protection checks pass.

Provenance, Python syntax, project contract and whitespace checks pass. Logs:
`build/eb-tree-selected-checks/{build,ctest,selected-ctest}.log`. Artifacts:
`build/selected-tree-mpi-debug/tests/selected_eb_tree_3d/`. Selected methanol,
selected adaptive schedules, MPI, installed binaries, convergence and physical
validation were not rerun in this increment.
