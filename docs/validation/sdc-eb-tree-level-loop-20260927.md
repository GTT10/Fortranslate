# Coupled EB tree SDC level update — 2026-09-27

## Source and behavior

Implementation `4b70f19a`, based on `088485a2` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The coupled serial EB SDC driver now accepts `target_depth`, ancestor endpoint
trees with start/end fractions, and optional ancestor source endpoint fields.
Hydro uses starting gas fractions and midpoint source fractions; old/new
transport stages use starting/ending gas fractions. Reaction and EOS recovery
operate only on the selected level. Other levels retain gas, temperature and
reaction estimates. Initial restriction and final reflux/restriction are
skipped for level updates, leaving synchronization to the recursive driver.

Level updates require an integrated-flux output. Optional old/new predictor
fields combine the corresponding external and transport rates with reaction
estimates. Inactive flux arrays remain unallocated and inactive predictors are
zero. Exchange inventories count only exposed cells on the selected level;
only a root-level update reports physical boundary flux. All-level calls keep
the existing common-step path. Predictor finiteness is checked before output
publication.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- `test_eb_tree_full_3d` built successfully.
- CTest `unit_eb_tree_full_3d`: 1/1 passed, 59.17 s test / 59.22 s total.
- New check: depth 2 of a three-level/five-patch tree, chemistry, all molecular
  transport flags and forcing enabled, two SDC iterations. Ancestor end states
  differ by a factor of 1.01; both interval fractions are 0.25. A reference
  explicitly constructs the same interpolated ancestor state.
  Selected gas, temperature, reaction estimates, x/y/z fluxes, old/new
  predictors and exchange diagnostics match exactly. Inactive gas/temperature/
  reaction fields remain unchanged; fluxes are absent and predictors zero.
- Existing common-loop conservation and state/diagnostic rollback checks remain
  in the same executable. New flux/predictor outputs were compared on success;
  their failure publication boundary was inspected, not separately fault-injected.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_full_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-eb-tree-level-loop-checks/`.
A recursive EB SDC trajectory, varying start/end ancestor fractions across a
coupled interval, selected build, MPI run, full matrix, installed-binary audit
and physical validation were not performed. Recursive orchestration remains
open; this check establishes its level-update interface only.
