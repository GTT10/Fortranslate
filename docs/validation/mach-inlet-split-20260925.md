# Mach inlet reactive composition — 2026-09-25

Implementation `45fc965`, based on `34675c0`.
`advance_reactive_full_3d` appends optional `inflow_mach` and
`inflow_pressure`, forwarding them through both molecular half steps and the
physical hydro step. Chemistry remains the existing outer pair of half steps.
Each boundary operator rebuilds its inlet from its current stage state.
Invalid active Mach controls on non-inflow faces are rejected before dispatch,
including periodic hydro with transport disabled. Pressure controls require a
Mach array and positive pressure requires an active Mach face. The existing
private candidate preserves fluid state and temperature on failure; diagnostic
outputs retain their existing intent(out) contract.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed unit target and selected H2O2 spray runtime
built successfully; selected execution was not run.
`unit_reactive_boundary_transport_3d`: **1/1 passed**, 8.68 seconds, total 8.72;
no failures or skips. The new composition check uses active elementary chemistry,
all molecular transport terms, characteristic PLM/HLLC, Mach 0.2 and prescribed
inlet pressure 120000 Pa. Conserved state, temperature and boundary inventory
match explicit R/2–T/2–H–T/2–R/2 calls exactly. Invalid inlet temperature after
the first chemistry stage leaves the caller's fluid unchanged; an active Mach
face cannot silently enter periodic dispatch. Existing six-face molecular and
wall checks also passed. No existing tolerances changed.

Mechanism provenance, Python syntax of the contract checker, project contract
and whitespace checks passed. This is serial composition evidence, not input-file,
LES, AMR or MPI integration, external nozzle parity or physical validation.

## Commands and artifacts

Build and test logs: `build/mach-inlet-split-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_boundary_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_boundary_transport_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The MPI-related environment settings above are host-local workarounds, not
application defaults. A final indentation-only edit followed the successful
build/test; it does not change Fortran tokens.
