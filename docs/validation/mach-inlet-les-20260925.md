# Mach inlet regular LES operators — 2026-09-25

Implementation `6bc4c51`, based on `d828fd0`.
`compute_les_rhs_3d`, `advance_les_ssprk2_3d`,
`compute_dynamic_les_boundary_rhs_3d` and `advance_dynamic_les_ssprk2_3d`
append optional `inflow_mach(6)` and `inflow_pressure(6)`. The static models
(Smagorinsky, WALE, Vreman) use the shared two-layer Mach ghost constructor;
Dynamic uses its filter-radius-plus-five donor width. Both SSPRK2 stages rebuild
inlet state from the current interior composition/pressure. Existing wall
transforms, pressure outlets and surface-flux accounting are retained.

RHS entry points reject invalid Mach/pressure controls before periodic dispatch.
Pressure controls require a Mach array; positive pressure requires an active
Mach inlet. Negative Mach retains ordinary reservoir behavior. State and
temperature commit only after successful integration. Dynamic's RHS, diffusion
and boundary diagnostics retain atomic publication; static LES retains its
existing intent(out) diagnostics rather than promising diagnostic rollback.
No model coefficients, numerical tolerances or timestep bounds changed.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed test/tree runtime and selected H2O2 spray
runtime built. Selected execution was not run.
`spray_les_dynamic_boundaries`: **1/1 passed**, 5.72 seconds, total 5.76.
The added checks cover all four models: a uniform initial Mach inlet RHS,
maximum diffusivity and surface rate exactly equal an independently constructed
reservoir oracle; nonuniform gas evolves under SSPRK2, with inventory change
matching its boundary budget at the existing 2e-14 scaled bound. Static cases
use extrapolated x-low pressure; Dynamic uses prescribed z-high pressure.
Invalid temperature preserves fluid state, Dynamic also preserves its boundary
diagnostic, and active Mach on a periodic face is rejected. Existing physical
Dynamic/AMR-root parity, wall and pressure tests pass in the same executable.

Initial run failed the new exact oracle check (4.55 seconds, total 4.60): the
fixture replaced the kernel's EOS-returned inlet temperature with literal
1100 K. Retaining the returned temperature fixed the oracle; no production
change or tolerance relaxation was needed. The final run had no failures/skips.
Mechanism provenance, Python syntax of the contract checker, project contract
and whitespace passed. This is regular LES operator evidence, not spray
application, AMR/MPI Mach integration, external nozzle parity or physical
validation.

## Commands and artifacts

Logs: `build/mach-inlet-les-checks/`, including initial/final CTest logs and a
copy of the final `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_dynamic_boundaries --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_dynamic_boundaries$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

MPI-related environment overrides are host-local workarounds, not defaults.
