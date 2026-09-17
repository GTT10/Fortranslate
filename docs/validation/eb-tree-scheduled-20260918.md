# Scheduled multilevel EB topology changes — 2026-09-18

Implementation: `558d12b438def2095637c160629f175aaa7e2558`, based on
`4366433`. Application: `pelef_reactive_eb_tree_3d`.

The namelist accepts up to 16 prescribed event times and complete replacement
topologies. Startup validates parent ordering, boxes, non-overlap and common
refinement ratios, then builds immutable geometry on the deepest scheduled
lattice. Each accepted step is clipped to the next event. Conservative remap
runs before budget validation and checkpoint saving. Output masking follows
the active topology rather than the initial input arrays.

V3 application envelopes bind initial and all scheduled topologies, event
times and completed stage. Restart initializes the expected active topology
before reading tree records and verifies the stage against saved time. Static
V1/V2 parsing remains, with schedules rejected if absent from the envelope.
No temperature tagging, recursive subcycling, MPI or moving geometry is added.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`; no toolchain configuration change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_eb_tree_application_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_eb_tree_(application|scheduled)_3d$' --output-on-failure --output-log build/eb-tree-schedule-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, zero failures/skips, total 9.20 s (static 2.55 s, scheduled
6.61 s). The schedule changes a three-level five-node tree to root only at
1.5e-9, then a three-level four-node tree at 2.5e-9. Four accepted steps reach
3e-9; output has **216 exposed fluid cells**. Maximum normalized inventory
residual is **7.05347627465518679e-15**, below the unchanged 5e-10 gate.
Independent CSV fluid-volume sums match final budgets.

Checkpoints after step one (before the first event), step two (root-only),
and step three (after refinement) all continue to byte-identical final state
and budget CSVs. Altered future event time and inconsistent stored stage reject
before output creation. Existing policy/mode and corrupted-record rejection
checks preserve sentinel output/checkpoint files. The static case also passes.
Automatic scheduled intervals, selected mechanisms, large-tree scaling and
physical/convergence qualification were not tested.

Provenance, Python syntax, project contract and whitespace checks pass. No
failed build/test attempt occurred. Logs:
`build/eb-tree-schedule-checks/{build,ctest}.log`. Runtime artifacts:
`build/selected-tree-mpi-debug/tests/eb_tree_scheduled_3d/`.
