# Regular 3D pressure outlet transport — 2026-09-18

Implementation `3301435`, base `c4045da`.

Hydro and physical transport now share `fill_pressure_outlet_ghosts_3d`.
Both normal ghost layers use the existing characteristic pressure kernel;
transverse ghost extents are included for transport gradients. Face order
defines corner composition. Optional `outlet_pressure(6)` is forwarded through
transport RK2 and from `advance_reactive_full_3d` to hydro and both transport
halves. Positive entries require an `outflow` face; zero retains extrapolation.
Invalid ghost states reject the private candidate before publication.

## Completed focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, Debug/Ninja, existing fixed tests-enabled
and selected H2O2 tests-disabled configurations. No MPI runtime test in this increment.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_boundaries_3d test_reactive_boundary_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_boundar(ies|y_transport)_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**, 0 skips, 11.11 s. Both builds and static checks
passed. The initial fixed build failed on overlong declarations; splitting
the declarations resolved it. No tolerance was changed.

The added transport fixture checks each face separately at target 102325 Pa
from uniform 101325 Pa, with viscosity, conduction, diffusion and barodiffusion
enabled. Nonzero exchange is required; integrated state change matches the
surface ledger at relative 3e-12 with scale floor 1e-20. A target producing
invalid ghost pressure rejects while preserving state and temperature.
Existing hydro and wall/transport regressions also pass.

Full reactive pressure forwarding was compiled, not independently exercised
with chemistry. Multiple pressure faces at corners, application input, AMR,
MPI, acoustic reflection and physical validation remain outside this evidence.
This is not a GC-NSCBC implementation.

Artifacts: `build/pressure-outlet-transport-checks/` contains `provenance.log`,
`build.log`, `build-final.log`, `build-selected.log`, and `ctest.log`.
