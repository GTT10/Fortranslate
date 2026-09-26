# EB tree SDC level stages — 2026-09-27

## Source and behavior

Implementation `f076b883`, based on `fbda2ba8` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Both provisional EB tree hydro and transport stages now accept zero-based
`target_depth` and ancestor start/end trees with interpolation fractions.
Hydro additionally accepts independent source start/end fields and fractions,
so gas and predictor time centering can be chosen separately. Ancestor inputs
require level selection. All-level calls retain their existing behavior.

Only selected-level patches are evaluated. Returned inactive rates are zero
and inactive flux arrays are unallocated. Existing transactional publication
is retained. The stages still perform no restriction, redistribution or reflux.
The coupled level loop and recursive orchestration are separate remaining work.

## Independently completed checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Provenance preflight passed for `mechanisms/h2o2_full.json`.
- Hydro and transport test targets built successfully.
- CTest: 2/2 passed, 48.83 s total (hydro 44.35 s, transport 4.43 s).
- New checks select depth 2 in a three-level/five-patch tree. Ancestor end
  states are scaled by 1.01; interpolation fraction is 0.25. A separate
  reference explicitly constructs the interpolated ancestor state before
  calling the stage. Selected rates and x/y/z integrated fluxes match exactly
  for both hydro and transport. Inactive rates are zero and fluxes absent.
- Hydro also forwards source endpoint/fraction arguments using identical
  nonzero source endpoints at fraction 0.5. Nontrivial source interpolation
  arithmetic is covered by the earlier source-halo checks, not newly varied
  in this stage regression.
- Existing common-step, conservation and rollback checks remain in both
  executables. No tolerances were relaxed.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_hydro_3d test_eb_tree_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(hydro|transport)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-eb-tree-level-stage-checks/`.
No recursive SDC trajectory, application restart, selected build, MPI run,
full matrix, installed-binary audit or physical validation was performed.
