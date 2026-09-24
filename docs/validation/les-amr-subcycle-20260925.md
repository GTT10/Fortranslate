# Recursive serial LES subcycling — 2026-09-25

Implementation `77f81eb`, based on `9abaad2`. The new
`advance_amr_patch_tree_subcycle_les_3d` driver uses the level-local LES API to
advance children with refinement-ratio-squared intervals. Coarse start/end
states supply temporal donor interpolation. Child interval fluxes accumulate
before conservative synchronization/reflux and temperature recovery. The whole
root interval publishes atomically, including boundary change, maximum reflux
correction, maximum diffusivity and optional per-level step counts.

The level-local final temperature recovery now uses the first-stage temperature
guess, matching the existing common-step LES path. No flux formula or tolerance
changed. The common-step API remains available. Automatic coarse-interval
selection and MPI/coupled application dispatch are not yet connected; each
actual level step still checks its diffusion bound and rejects an excessive dt.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, existing fixed tests-enabled and selected H2O2
tests-disabled configurations. Fixed test/runtime and selected runtime built.
Selected execution and MPI subcycling were not run.

- Final `spray_les_amr_subcycle`: **1/1 passed**, 83.45 seconds, total 83.49.
  Smagorinsky, WALE, Vreman and Dynamic closures exercise root-only and nested
  trees, including peers and refinement ratios 2 and 3. Counts are exactly
  1/4/36; composite conservation is within `2e-10` relative to inventory,
  parents equal independent restriction exactly, temperatures remain positive,
  and species closure is within `1e-12`. Periodic, moving/isothermal-wall and
  pressure-outlet budgets pass. Excessive dt preserves gas, temperatures,
  budget, correction, diffusivity and step-count outputs.
- `spray_les_amr_level` after the recovery-guess change: **passed**, 2.75 seconds.
  Root common-step parity, temporal coarse donors, integrated flux accounting,
  inactive-level preservation and rollback remain covered.

Two development failures were resolved. The first recursive test failed the
root temperature comparison (0.21 seconds): the level API used an averaged
temperature guess, unlike the common-step path. Aligning that guess fixed it
without changing the `1e-9` K recursive-test tolerance. The next run passed
the level test but rejected the Dynamic root pressure fixture (recursive test
13.55 seconds): two-direction seven-layer acoustic corner extrapolation made
that fixture inadmissible. The Dynamic fixture now prescribes pressure on one
outlet axis; the three static fixtures retain two prescribed outlet axes.
Pressure admissibility checks remain unchanged. The final recursive run passed.

No skips occurred. Provenance preflight, Python syntax of the contract checker,
project contract and whitespace checks passed. These are focused integration
checks, not external PeleC integrator parity, convergence, physical validation
or a full configuration matrix.

## Commands and artifacts

Build logs and initial, intermediate and corrected CTest logs are in
`build/les-amr-subcycle-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_subcycle test_les_amr_level --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_amr_(subcycle|level)$' --output-on-failure
# After correcting the Dynamic corner fixture, rerun the affected test:
ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_amr_subcycle$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
