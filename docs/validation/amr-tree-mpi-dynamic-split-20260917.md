# Dynamic recursive MPI full-split integration — 2026-09-17

## Source and purpose

Base `3a7d33e08907eb8593b3d66ff041834e35773dd9`, local `main`.
Adds a focused application regression for the combined scheduled-regrid,
recursive reaction/transport/hydro, measured-work rebalance and indexed restart
path. Solver code, mechanisms, tolerances and golden hashes are unchanged.

`tools/check_mpi_amr_tree_dynamic_split_3d.py` SHA-256:
`59f283642a332afd34f67fc22beecb2b73ede6b8c217c5a06b4bb3f370711c9a`.

## Fixture

The existing 4-cubed elementary hotspot input retains chemistry, viscosity,
conduction, species diffusion and prescribed tangential wall motion. It selects
characteristic PLM/MC with recursive stepping. Small generated topologies contain
four nodes/three levels initially, three nodes/three levels after a move and
ratio change (2 to 3 at the first refined depth), then only the root, and finally
the initial branched tree again. Events at 2.5e-10, 5e-10 and 7.5e-10 s precede
final time 1e-9 s. All three events are followed by an integration interval.

For each fixed/selected elementary frontend:

1. Serial continuous execution supplies final CSV/budget bytes over four steps.
2. Two ranks with weighted ownership and every-step measured-work rebalancing
   save indexed V2 at step one, after the first geometry change.
3. Three ranks with round-robin restart ownership and rebalancing continue across
   coarsening to the root (two empty ranks), refinement back to three levels and
   a final integration interval; final CSV/budget must equal serial bytes.
4. Serial step-one save must equal the MPI checkpoint bytes exactly and must not
   overwrite the already completed output files.

Fixed/selected continuous bytes must also match. This is eight small application
runs in one CTest, not eight independent CTests. Timing-driven owner choices are
not prescribed, and no speedup assertion is made.

## Commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6; existing matching build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d -j 4
python3 -m py_compile tools/check_mpi_amr_tree_dynamic_split_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_dynamic_split_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/amr-tree-mpi-dynamic-split-checks/{build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_dynamic_split_3d/`,
including generated geometry/schedule, input files, run logs, final CSV/budget and
V2 checkpoints for each mechanism frontend.

## Completed results

Build and CTest **1/1 passed**, no failures/skips, 106.05 s total (test 106.02 s).
Provenance, Python syntax, project contract and whitespace checks passed.
Both mechanisms completed four steps and three scheduled regrids. Maximum reported
remap error was 5.2943e-18, maximum reflux 6.8795e-5 and minimum transport theta 1.
Both three-rank continuations actually migrated owners at step four. Serial/MPI
final CSV, budget and step-one V2 checkpoint byte comparisons all passed, as did
fixed/selected equality. No solver changes were necessary.

## Limits

Prescribed geometry changes only; full-split temperature-triggered adaptation,
large trees, production scaling, detailed-mechanism CFD, installed audits,
external PeleC parity and physical validation remain unqualified. The short small
case does not establish time/space convergence or a physically calibrated flame.
