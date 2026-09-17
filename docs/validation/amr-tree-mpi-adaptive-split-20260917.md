# Temperature-adapted recursive MPI full split — 2026-09-17

## Source and scope

Base `596fdd56e3245d0389cca8b64960dc53b4eba105`, local `main`.
Adds `--adaptive` to the existing combined dynamic-split regression. This closes
the focused integration check for temperature-driven geometry, recursive
reaction/transport/hydro, per-step measured-work rebalancing and indexed V2
changed-rank continuation. Solver formulas, mechanisms, tolerances and golden
hashes are unchanged.

`tools/check_mpi_amr_tree_dynamic_split_3d.py` SHA-256:
`5ce4211d9029c8cc5ae986d3c52ae132926693948514b08ca4e712bf37cb5253`.

## Fixture and assertions

The 4-cubed hotspot starts with only the root and an off-center 1200 K peak.
Chemistry, viscosity, conduction, species diffusion and tangentially moving walls
remain enabled; characteristic PLM/MC and recursive stepping are selected.
Every completed step checks temperature thresholds 1150 K on two refinement
levels with ratios 2/2, single-cell tiles and no tag buffer. CFL is 0.001 and
final time 1e-9 s. The maximum node budget is 64.

For each fixed/selected elementary mechanism, serial continuous output supplies
the reference, a weighted two-rank run saves after its first adaptation, and
three ranks resume with round-robin restart owners and per-step work rebalancing.
The checker verifies step/adaptation counts, no repeated first-step tag on restart,
final level 2, exposed-cell volume, exact final CSV/budget bytes and exact serial/MPI
V2 checkpoint bytes. It also compares fixed/selected bytes. Serial stop must leave
completed final outputs intact. Eight small application runs constitute one CTest.

## Commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing matching build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d -j 4
python3 -m py_compile tools/check_mpi_amr_tree_dynamic_split_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_adaptive_split_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/amr-tree-mpi-adaptive-split-checks/{build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_adaptive_split_3d/`,
including each mechanism's inputs, logs, final CSV/budget and indexed checkpoints.

## Completed results

Build target check passed (`ninja: no work to do` after CMake regeneration).
CTest **1/1 passed**, no failures/skips, 162.59 s total (test 162.57 s).
Provenance, Python syntax, project contract and whitespace checks passed.
Both mechanisms refined from one root to ten nodes/three levels at step one,
then completed three steps with three tag checks. The maximum reported remap
error was 2.1071e-16. Both three-rank continuations migrated owners at steps two
and three. Final CSV/budget and step-one V2 checkpoint comparisons were byte-exact
against serial; fixed/selected bytes also matched. No solver changes were needed.

## Limits

The hotspot triggers refinement; it does not subsequently force temperature-driven
coarsening or patch motion. Those geometry transitions have separate prescribed-
schedule evidence. Timing-driven placements are not fixed expectations. Detailed
mechanism CFD, large-tree performance, convergence, full configuration matrices,
installed audits, external PeleC parity and physical validation were not run.
The scheduled branch of this checker was qualified in the preceding commit and
was not independently rerun here.
