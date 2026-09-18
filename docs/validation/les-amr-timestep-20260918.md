# Serial AMR LES timestep — 2026-09-18

Implementation `e698432e7943ac4fc737bfacc1e04a7f73ec0a06`, based on `763c610`.

`compute_amr_patch_tree_les_timestep_3d` selects the minimum of
`cfl/(D_patch*sum(1/h_patch^2))`, capped at `1e6` seconds. CFL must be finite
and in `(0,0.2]`. It averages covered cells and recovers temperature on private
copies, packs the tree, and uses the same coarse/peer/physical ghost filling,
pressure outlet correction and static LES RHS as the first integration stage.
Moving/isothermal walls, optional coarse limiting and CI enter that evaluation.
The maximum returned diffusivity is diagnostic; grid spacing determines the
limiting patch. No caller state is modified, and timestep/diffusivity publish
only after all patches succeed. Zero diffusion returns the cap.

This is an initial-state frozen-gradient estimate. Both integration stages
retain their own interval and EOS checks; callers must handle rejection and
combine this bound with other operators and event times. MPI, subcycling and
application use of this selector remain pending.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug, unchanged configurations.
Fixed tests-enabled `test_les_amr_tree` and selected H2O2 tests-disabled
`pelef_selected_spray_tree_runtime` built. Selected runtime was not run.
CTest **1/1 passed**, 4.52 seconds (4.56 s total), zero failures/skips.

The existing root/six-node AMR test now exercises the selector for all three
static closures with CI, periodic/outflow/moving-isothermal/pressure boundaries
and MC coarse interpolation. Root CFL agrees with direct regular RHS diffusion
within `1e-12`; halving CFL exactly halves dt and leaves diffusivity unchanged.
Invalid CFL 0.21 preserves output sentinels. Steps use the minimum of the selected
interval and `1e-9 s`; this does not demonstrate acceptance at the full bound.
The preceding root state/temperature/ledger parity, composite conservation and
excessive-interval rollback assertions also pass, without tolerance changes.
Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, full matrix, MPI execution, external solver parity or
physical validation occurred in this increment.

## Commands and artifacts

Logs: `build/les-amr-timestep-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_tree --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_amr_tree$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
