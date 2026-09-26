# Composite redistribution in common-step EB SDC — 2026-09-27

Implementation `cc1ec227`, based on `d1792f3a` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`advance_reactive_eb_tree_sdc_3d` now accepts optional
`redistribution='flux_redist'`. The default remains `none`. Common-step mode
redistributes old/new molecular-transport rates and the hydro rate over the
composite leaf mesh. The same linear redistribution is applied to the final
reflux increment, followed by temperature recovery and restriction. Reaction
and forcing inventories retain their existing definitions. Fluxes still carry
the conservative transfers used for synchronization and boundary accounting.

The source-rate operator additionally accepts boundary labels. Paired periodic
faces wrap integer neighbor coordinates, including when opposite faces have
different refinement levels. Boundary-label and root periodic-geometry validation
precede output publication. Other domain faces remain closed for redistribution.

Selected-level redistribution is explicitly rejected: a composite neighborhood
needs rates at a common stage time. Recursive transfer accounting and MPI
redistribution remain future work. Application input does not yet expose this
library selector, and state redistribution is not implemented for the tree.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`; existing compiler warnings remain.

- Mechanism provenance and both touched test targets built successfully.
- CTest: 2/2 passed, 13.44 s total; operator test 0.02 s, coupled SDC test 13.38 s.
- Operator checks now include translation by one cell through both a periodic
  domain seam and a peer-patch seam, with absolute discrepancy below 2e-14.
  Coarse-to-fine and fine-to-coarse periodic impulses reach the opposite side
  and conserve signed components within 3e-14 absolute integrated discrepancy.
  Unpaired periodic labels and mismatched periodic geometry reject. Output
  retention is checked after an unpaired-label failure.
- Coupled test: three levels, five patches, all transport components, reactions,
  external force, two iterations and dt 1e-10 s. The redistributed state differs
  from the baseline; relative inventory residual is 4.0936429451177748e-15
  (`none`: 3.7650552248896492e-15). Chemistry and forcing remain active.
- Invalid artificial viscosity after initial transport evaluation rejects while
  preserving gas, temperature, reaction estimates, inventories and limiter.
- Python syntax, project contract and `git diff --check` passed. No final build
  or test failures, skips or interrupted runs were counted.

Artifacts: `build/sdc-eb-tree-redist-periodic-checks/`, with `provenance.log`,
`build.log`, `ctest.log`, `detail.log` and `contract.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_eb_tree_full_3d test_eb_tree_redistribution_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_eb_tree_(sdc_)?redistribution_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

This is bounded operator and library-stage evidence. The final flux synchronization
still checks the intermediate unredistributed reflux state for admissibility;
redistribution does not rescue a rejected intermediate state. No larger stable
time step, temporal order, upstream AMR redistribution parity, selected-mechanism
build, installed-binary audit or physical validation is established here.
