# EB tree SDC regrid/restart — 2026-09-27

## Source and behavior

Implementation `7ece3707`, based on `67521afc` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`remap_reactive_eb_tree_3d` optionally remaps the packed species/energy reaction
estimate together with gas. It first restricts source descendants with fluid
volume weights, maps same-level or nearest ancestor PCM donors, then restricts
the new hierarchy. Signed values are not clipped or interpreted as gas. Solid
target cells receive zero. Both gas and source publish only after the entire
remap succeeds. The application forwards its persistent SDC source and refreshes
the checkpoint record after activation of the new topology.

The previous scheduled/adaptive SDC regrid gate is removed. Common timesteps
remain required; MPI SDC remains gated.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Provenance preflight passed for `mechanisms/h2o2_full.json`.
- Application and `test_eb_tree_remap_3d` targets built successfully.
- CTest `unit_eb_tree_remap_3d`: 1/1 passed, 1.02 s test / 1.07 s total.
  Added signed source movement, complete coarsening and refinement checks;
  composite source integrals satisfy 3e-14 times max(1, abs(initial)).
  Invalid gas causes both gas and source outputs to remain unchanged.
- Coupled application: three-level/five-patch initial topology, coarsen to
  one root at 1e-10 s, refine to five patches at 2e-10 s, finish at 3e-10 s.
  Chemistry, transport and forcing are enabled. Continuous three-step run
  versus stop/checkpoint after first regrid and restart: final state and budget
  CSVs are byte-identical. Both regrid stages are observed in the logs.
- Continuous composite inventory error: 8.3890849755363664e-15.
- Python syntax, project contract and `git diff --check` passed.

Artifacts and trajectory runner: `build/sdc-eb-tree-regrid-checks/`.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_remap_3d pelef_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_remap_3d$' --output-on-failure
python3 build/sdc-eb-tree-regrid-checks/run.py
python3 -m py_compile tools/check_project_contract.py build/sdc-eb-tree-regrid-checks/run.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The application trajectory exercises scheduled regridding. Temperature-based
adaptive tagging shares activation but was not independently run. No selected
build/run, MPI run, recursive SDC step, full matrix, installed-binary audit or
physical validation was performed. PCM remapping preserves the tested volume
integrals; this does not establish spatial convergence or full PeleC parity.
