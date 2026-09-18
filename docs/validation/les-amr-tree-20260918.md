# Serial AMR static LES — 2026-09-18

Implementation `090237b58434ab86a7317c0e585f4c70dd92131b`, based on `51386d5`.

`les_amr_patch_tree_mod::advance_amr_patch_tree_les_3d` connects the shared static
LES face kernel to the established AMR ghost/reflux APIs. The orchestration
follows the common-step molecular-transport integrator: average initial covered
cells, recover temperatures, fill peer/coarse/physical halos, evolve each patch,
synchronize/reflux accepted face transfers and recover EOS after each Euler
stage, then form the SSPRK2 combination and recover temperature. Flux integrals
are dt times flux density; existing reflux handles area weighting and composite
boundary corrections. A failed stage preserves caller state, temperature,
boundary ledger and maximum-correction diagnostic.

The final EOS recovery uses the first-stage temperature guess, as the regular
LES SSPRK2 path does. Supported closures are Smagorinsky, WALE and Vreman with
optional CI. Coarse donor limiting, moving/isothermal walls, reservoir inflow
arguments and pressure outlets use existing shared helpers. Positive dt is
required; each patch/stage checks the existing diffusion-number bound 0.2.
This is a frozen-gradient heuristic, not a nonlinear stability proof. MPI,
subcycling, automatic hierarchy dt selection, dynamic LES, EB and application
integration remain open. The optional spray tree runtime contains the module.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled `test_les_amr_tree` and selected H2O2 tests-disabled
`pelef_selected_spray_tree_runtime` built. Selected runtime was not run.
CTest **1/1 passed**, 2.50 seconds (2.54 s total), zero final failures/skips.

The test spans all three closures with CI 0.1, a root-only grid and a six-node
hierarchy containing neighboring patches and four nested refinement levels.
Each uses periodic, outflow and moving-isothermal-wall/pressure-outlet boundaries
with MC coarse interpolation and anisotropic root spacing. Root-only state
matches the regular LES API within `2e-13` relative with unit floor, temperature
within `1e-9 K`, boundary ledger within `1e-14`. Composite mass/momentum/energy/
species inventories match boundary changes within `2e-12` relative with `1e-12`
scale floor. Overlarge dt rejects with exact state/temperature/diagnostic rollback.

The first run failed root temperature parity when final recovery used the
molecular-transport integrator's averaged temperature guess. The final guess
was aligned with regular LES; the same test/tolerances then passed. Initial
failure is retained in `ctest-temperature.log`. Provenance, Python syntax,
project contract and whitespace checks passed. No full matrix, MPI execution,
external solver parity, convergence study or physical validation ran.

## Commands and artifacts

Logs: `build/les-amr-tree-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `ctest.log`, `ctest-temperature.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_tree --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_amr_tree$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
