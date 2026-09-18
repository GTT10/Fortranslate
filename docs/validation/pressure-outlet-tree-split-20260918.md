# Patch-tree pressure outlet subcycled transport and split update — 2026-09-18

Implementation `b9c6207`, base `65b3c2d`.

The serial recursive transport API forwards optional `outlet_pressure(6)` to
all level updates. The full R/2–T/2–H–T/2–R/2 wrapper forwards it to both
transport halves and hydro, with common or recursive timesteps. Existing
private candidates retain atomic publication of state and diagnostics.
MPI tree forwarding, standalone reaction/hydro Strang entry points and tree
application input remain pending.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, reused Debug/Ninja fixed tests-enabled and
selected H2O2 tests-disabled configurations; no MPI runtime test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_transport_subcycle_3d test_amr_patch_tree_split_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_(transport_subcycle|split)_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Transport cases add x-high/z-high targets of 102325 Pa with moving isothermal
y walls, root-only and four-node three-level trees. Checks retain root common
step comparisons, 1/4/36 recursive step counts, composite conservation,
species closure, synchronized parents and large-step transactional rejection.

Full split cases add 203650 Pa targets on x-high/z-high for common explicit
chemistry and recursive implicit chemistry. The four-node tree checks Euler
and element inventories against the boundary ledger, species closure,
restriction, post-reaction failure rollback and disabled chemistry/transport
hydro parity. Existing zero-pressure comparisons remain. Tolerances are unchanged;
these are implementation checks, not physical validation.

Artifacts: `build/pressure-outlet-tree-split-checks/{provenance,build,build-selected,ctest}.log`.

Final CTest **2/2 passed**, zero skips, 106.79 s total (split 83.74 s,
recursive transport 23.01 s). Fixed/selected builds, provenance, Python
syntax, project contract and diff checks passed. No failed or interrupted runs.
