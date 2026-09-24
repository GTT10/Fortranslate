# LES subcycle timestep selection — 2026-09-25

Implementation `1752e1f`, based on `82bc013`.
`compute_amr_patch_tree_subcycle_les_timestep_3d` selects a coarse interval for
the existing serial recursive LES driver. It reuses the common selector's
averaged state, ghost donors, physical boundaries and model diffusivity.
For level step count N, its rate is `D * sum(1/h**2) / N`; N is the cumulative
product of refinement ratios squared, checked for int64 overflow.
The smallest coarse-interval bound across all patches is returned.

The common selector gains optional `subcycle=.true.`; its default retains
N=1. Both selectors retain CFL range `(0,0.2]`, the `1e6` interval cap,
maximum diffusivity diagnostic and transactional output publication.
The recursive driver still checks every stage against its evolving diffusion
bound. Selection from the initial state does not guarantee that a later
nonlinear stage will be admissible.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed tests-enabled runtime/test and selected
H2O2 tests-disabled runtime built. Selected execution was not performed.
Focused CTest `spray_les_amr_level`: **1/1 passed**, 4.98 seconds, total 5.02.

All four LES models verify root-only equality with common selection, nested
interval bounds between common dt and four times common dt, unchanged maximum
diffusivity, exact proportional response to halving CFL, and output preservation
for invalid CFL. The Dynamic fixture also verifies a strictly larger coarse dt
when its child-level restriction is scaled. Existing level-step parity,
temporal coarse donor, integrated-flux and rollback assertions remain covered.

No failures, skips or tolerance changes occurred. Mechanism provenance,
Python syntax of the contract checker, project contract and whitespace passed.
This is selector and level-operator evidence. The recursive conservation test
was not repeated; MPI selection, MPI recursion and application dispatch remain
implementation work. No external solver or physical validation is claimed.

## Commands and artifacts

Logs: `build/les-subcycle-timestep-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_level --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_amr_level$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
