# Mach inlet molecular transport — 2026-09-25

Implementation `5271816`, based on `d06e6ee`.
`physical_transport_fluxes_3d`, `reactive_transport_euler_update_3d` and
`advance_reactive_transport_3d` append optional `inflow_mach` and
`inflow_pressure`. Physical flux evaluation uses the shared Mach inlet ghost
constructor with viscous wall transforms. Both SSPRK2 stages rebuild inlet
pressure/composition from that stage's state, retaining molecular viscosity,
conduction, diffusion, barodiffusion, donor limiting and boundary accounting.
The absent-control path retains its existing reservoir behavior.

The Euler entry rejects active Mach controls on non-inflow faces and invalid
Mach/pressure values instead of silently taking a periodic path. Positive inlet
pressure requires an active Mach face; pressure controls require a Mach array.
Mach-only inflow does not need a fixed reservoir-state table. The existing
advance routine still returns early for zero intervals or disabled transport;
this change does not redefine that no-op input-validation contract. State and
temperature still publish only after both stages succeed; existing intent(out)
diagnostic behavior is unchanged. Input-file, splitting, LES, AMR and MPI
forwarding remain to be connected.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed core/unit and selected H2O2 spray runtime
built. Selected execution was not run. `unit_reactive_boundary_transport_3d`:
**1/1 passed**, 7.39 seconds, total 7.43. No failures or skips.

All six faces exercise extrapolated or prescribed pressure. For a uniform
initial interior, a Mach Euler step equals its explicitly constructed reservoir
oracle exactly in conserved state, temperature, limiter and boundary inventory.
A full SSPRK2 step evolves the gas, preserves nonnegative species and species
closure, and matches inventory change to boundary flux at the existing 2e-12
scaled budget bound. All molecular terms, including barodiffusion, are enabled.
Invalid inlet temperature preserves gas and temperature. Existing wall,
pressure-outlet and molecular boundary checks also pass in the same test.
No existing tolerances were changed. Mechanism provenance, Python syntax of the
contract checker, project contract and whitespace passed. This is molecular
boundary-operator evidence, not application/MPI execution, external nozzle parity
or physical validation.

## Commands and artifacts

Logs: `build/mach-inlet-transport-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_boundary_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_boundary_transport_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
