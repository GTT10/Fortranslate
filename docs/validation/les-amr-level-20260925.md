# Level-local LES integration — 2026-09-25

Implementation `0fb4868`, based on `af95f25`. The serial API
`advance_amr_patch_tree_les_level_3d` advances only the selected tree depth with
SSPRK2, using the existing Smagorinsky, WALE, Vreman or Dynamic flux kernel.
Coarser packed start/end states and per-donor-level stage fractions supply
time-interpolated ghosts. Endpoints contain conserved components followed by
one temperature guess, following the existing transport-level contract.

Successful publication returns state/temperature, the maximum diffusion bound,
and half-dt-weighted sums of the two stage face fluxes. Only active-depth flux
entries are allocated. Inactive state and temperature remain unchanged. All
outputs remain unchanged on failure, including pre-existing integrated fluxes.
Each stage enforces the existing LES diffusion interval bound.

The recursive caller must accumulate child-interval fluxes and perform
restriction/reflux. This API deliberately does not synchronize time levels.
The common-step API is unchanged. Recursive scheduling, MPI level-local LES
and application dispatch remain required; this is not completed subcycling.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed tests-enabled runtime/test and selected
H2O2 tests-disabled runtime built. Selected execution was not performed.

Final CTest `spray_les_amr_level`: **1/1 passed**, 2.69 seconds, total 2.73.
All four models exercise:

- Exact root-only state agreement with the existing common-step operator,
  and temperature agreement within `1e-8` K.
- Active child evolution with changing coarse endpoints and MC reconstruction;
  a frozen-endpoint control produces a different accepted child state.
- Exact preservation of inactive fields and allocation only of active fluxes.
- Accepted state change accounted for by integrated face divergence, within
  `2e-12` times the initial state maximum magnitude.
- Missing endpoints and excessive dt rejection preserving state, temperature,
  diffusion output and all three integrated flux arrays.

An earlier test revision passed in 1.68 seconds before adding the changing-
endpoint control and stage-failure checks. No failures, skips or tolerance
relaxations occurred. Mechanism provenance, Python syntax, project contract
and whitespace checks passed. No MPI, full-matrix, external solver parity or
physical validation is claimed.

## Commands and artifacts

Build logs and final `LastTest.log` are in `build/les-amr-level-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_level --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_amr_level$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
