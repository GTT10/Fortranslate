# Recursive tree molecular-transport subcycling — 2026-09-17

## Source and behavior

Base: `bdde431c05dbe29d8a4b08382c1c0d72e08f01d2`, local `main`.
The transport module gains a per-depth SSPRK2 primitive returning integrated
fluxes and limiter diagnostics. A new scheduler advances children with ratio-squared
substeps, interpolates coarse endpoints, reconciles interval fluxes, recovers
temperature and commits the complete coarse interval transactionally. The timestep
helper scales local diffusion and touched-isothermal-wall bounds by cumulative
ratio squared while retaining the actual spatial spacing.

SHA-256 at verification:

- `src/amr/amr_patch_tree_transport_3d_mod.F90`:
  `5ff1eb9928c54a61baa0128bb4350129ac0e69088cdadbb8cf3e3658051f304a`
- `src/amr/amr_patch_tree_transport_subcycle_3d_mod.F90`:
  `a7448aa2a1f27eabfd27c53ee135d891bc3e6d2e30ca50ba29ba89b4ed11a333`
- `tests/unit/test_amr_patch_tree_transport_subcycle_3d.F90`:
  `1d7e039ec42bb41e2cc3da4706f94c024f25e9876146227d42b3685ce2e1f9ed`

## Completed checks

GNU Fortran 13.3.0 Debug, bounds checks/backtraces/FPE traps; MPI disabled.
Existing build: `build/selected-amr-tree-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/selected-amr-tree-debug --target pelef_core pelef_selected_amr_tree_reactive_3d_runtime -j 6
cmake --build build/selected-amr-tree-debug --target test_amr_patch_tree_transport_subcycle_3d -j 6
ctest --test-dir build/selected-amr-tree-debug -R '^unit_amr_patch_tree_transport_subcycle_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed and selected runtime builds succeeded. Focused CTest: **1/1 passed**,
14.80 s, no failures or skips. Provenance, project contract and whitespace passed.
No Python repository source changed; Python syntax verification was not applicable.

The unit uses root-only and branched three-level trees with ratios 2 then 3,
all molecular transport terms active, MC coarse donors, periodic boundaries or
outflow plus moving isothermal walls. It checks 1/4/36 step counts, scaled
coarse diffusion bounds, root agreement with the existing synchronous transport
solver, composite conservation including boundary flux, exact parent restriction,
positive temperatures, species closure and limiter bounds. A rejected large step
preserves state, temperature, budget, maximum correction, theta and step counts.

Logs: `build/amr-tree-transport-subcycle-checks/`.
Executable: `build/selected-amr-tree-debug/tests/test_amr_patch_tree_transport_subcycle_3d`.

## Limits

These are serial transport APIs. Coupled reaction/transport/hydro dispatch and
application timestep/checkpoint controls remain to be connected; tree MPI remains
open. No convergence order, performance/scaling, compiler matrix, installed-binary
audit, external PeleC parity or physical validation is claimed. The stepper accepts
caller-chosen dt and does not enforce the helper's bound or retry internally.
Reservoir bounds are still the caller's responsibility. Full-tree scratch copies
remain, and component interpolation does not itself prove EOS admissibility.
