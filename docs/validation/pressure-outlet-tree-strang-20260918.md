# Patch-tree pressure outlet reaction/hydro Strang APIs — 2026-09-18

Implementation `ebb3c46`, base `15e1587`.

Serial and MPI standalone reaction/hydro Strang entry points now accept optional
`outlet_pressure(6)` and forward it to hydro between reaction half steps. This
completes pressure forwarding for these APIs in addition to the already connected
full reaction/transport/hydro application path. Existing private-candidate
publication and collective hydro pressure validation remain unchanged.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, `/usr/bin/mpiexec`; reused Debug/Ninja fixed
tests-enabled and selected H2O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_chemistry_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_pressure_strang_3d_2$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_tree_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

A --pressure-only mode in the owner-local chemistry/Strang test adds a distributed three-node
pressure case: x-low/z-high targets 203650 Pa, slip y walls, characteristic PLM
and implicit chemistry over 1e-10 s. Two ranks compare every owned state/temperature entry
and boundary/reflux diagnostic exactly against serial, including active leaf
chemistry. Post-reaction hydro failure preserves the entire coupled state and
diagnostics. Invalid reaction and interval disagreement checks also run in the pressure case.
The default invocation retains the original three zero-pressure cases; these
were not successfully requalified in this turn. No tolerances changed.

Initial fixed build rejected an overlong test call; its line was wrapped and the
build rerun successfully. An initial combined four-case run timed out at
120.23 s. The new pressure case was isolated under its own CTest registration
with a short reacting interval; original case intervals and timeouts remain
unchanged. Selected compilation is separate from fixed-mechanism
runtime checks. This is implementation parity, not physical validation.

Artifacts: `build/pressure-outlet-tree-strang-checks/{provenance,build,build-final,build-focused,build-selected,ctest,ctest-focused}.log`.

Final focused CTest **1/1 passed**, zero skips, 10.18 s total (test 10.14 s).
Fixed/selected builds, provenance, Python syntax, project contract and diff
checks passed. The earlier combined-test timeout is not counted as a pass.
