# EB patch-tree SDC reaction stage — 2026-09-27

## Source and scope

Implementation `6a18ade4`, based on `7e715f6d` on local `main`.
PeleC reference remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`react_sdc_eb_tree_3d` accepts an old EB tree and packed nonreacting rates,
returning packed conserved state, temperature and species/energy reaction
estimates. The optional zero-based `target_depth` advances one level while
retaining caller outputs at other levels. All-level mode initializes gas from
the old tree. Solid cells retain old gas and caller reaction estimates.

The stage calls the existing forced SDC reactor only in fluid cells. Refined
parent cells also evolve provisionally: no averaging or reflux occurs inside
this stage. State, temperature and reaction outputs publish together only after
all selected cells succeed. Source arrays use the existing generic AMR tree
field container, with shape validation that permits unused solid-cell values.

## Independently completed checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- Target `test_eb_tree_full_3d` built successfully.
- CTest `unit_eb_tree_full_3d`: 1/1 passed, 45.93 s test / 45.98 s total.
  The existing split-step tests remain in this executable.
- New checks use three levels, five patches and a z-plane cut at 0.31.
  Forced reaction updates match the single-patch EB SDC reaction exactly for
  state, temperature and reaction estimates, including refined parents.
  Each of the three selected levels matches the corresponding all-level
  result; other levels retain sentinel outputs exactly.
- NaN nonreacting rates in solid cells are ignored. A negative density in the
  last fluid cell of the last patch rejects the stage after earlier patches
  have been evaluated, preserving every caller output exactly.
- Python syntax, project contract and `git diff --check` passed.

Commands from the repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_full_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-eb-tree-reaction-checks/{provenance,build,ctest}.log`.

This is a library reaction-stage check. Coupled tree SDC hydro/transport,
coarse/fine source halos, application/restart routing and MPI integration are
not established by it. No full configuration matrix, selected build,
installed-binary audit or physical validation was run in this step.
