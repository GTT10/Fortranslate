# Mach inlet ghosts and serial hydro — 2026-09-25

Implementation `ce65523`, based on `8f7260b`.
`fill_mach_inlet_boundary_ghosts_3d` applies the nozzle state kernel to regular
boundary halos (width 1..128, default 2). A six-entry Mach array enables this
behavior at existing `inflow` faces with nonnegative values; negative entries
retain ordinary reservoir/other face behavior. Active faces require prescribed
inlet temperature. Optional pressure zero selects interior pressure; positive
pressure selects the target. A Mach-only inlet needs no reservoir state table.
Ordinary reservoir faces retain their existing input requirements.

Each active ghost takes pressure and composition from its mapped interior
donor, recomputes the prescribed-temperature/Mach state, then applies the shared
wall reflection and tangential-motion transformation for corners. Existing axis
order determines inlet precedence at corners. Ghost state and temperature
publish together after success. The wrapper does not yet forward supplied
velocity fluctuation/total-velocity modes from the lower-level kernel.

`advance_reactive_boundary_hydro_3d` appends optional `inflow_mach` and
`inflow_pressure`. Both SSPRK2 stages rebuild these ghosts from that stage's
state before PCM or characteristic PLM fluxes. The existing surface ledger
accounts for inlet exchange. State/temperature publication remains transactional;
existing intent(out) diagnostic-zeroing behavior on failure is unchanged.
Absent Mach controls retain the original reservoir path. Input files, transport,
LES, AMR and MPI dispatch still require subsequent integration.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed core/unit and selected H2O2 spray runtime
built. Selected execution was not run. Final `unit_reactive_boundaries_3d`:
**1/1 passed**, 7.79 seconds, total 7.83, no skips.

Six faces verify exact seven-layer values against the state kernel and unchanged
interior cells. PCM with extrapolated pressure and characteristic PLM with fixed
pressure both produce positive net mass inflow; all conserved inventory changes
match surface fluxes within the existing 3e-12 scaled bound. Negative inlet
temperature rejects without changing gas or temperature. Previous kernel and
reservoir/pressure-boundary checks remain in the same test.

The initial fixed test build failed because two new integral calls omitted their
three grid-size arguments. The test calls were corrected; the final build and
execution passed. No numerical implementation change or tolerance adjustment was
needed. Mechanism provenance, Python syntax of the contract checker, project
contract and whitespace passed. This is serial boundary/hydro evidence, not
application execution, a frozen C++ nozzle comparison or physical validation.

## Commands and artifacts

Logs: `build/mach-inlet-hydro-checks/` (initial failed fixed build, final fixed
build, selected build and CTest logs).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_boundaries_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_boundaries_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
