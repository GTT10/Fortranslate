# Dynamic LES wall heat-flux kernel — 2026-09-25

Implementation `3e200fb`, based on `e736a30`. The ghosted Dynamic LES API accepts optional patch-local
`boundary_faces` and `wall_temperature`. The caller still supplies valid
velocity/state donor halos. Only patch faces identified as walls receive a
thermal correction; coarse/fine and interior faces must remain non-wall.

Adiabatic slip/no-slip walls set normal SGS heat flux to zero. Isothermal walls
use the face Dynamic thermal coefficient, arithmetic face strain viscosity and
mixture heat capacity with the half-cell temperature gradient, matching the
existing static LES boundary discretization. Stress and its face-velocity work
are preserved. The thermal diffusion bound includes the half-cell factor of
two. Wall corrections enter the returned faces and their conservative divergence
before transactional publication. Periodic calls need no new arguments.

The full padded-workspace routine also accepts explicit interior `wall_bounds`;
these are supplied automatically by the ghosted API. This is a local kernel
extension. AMR steppers and applications still reject Dynamic wall faces until
their dispatch and boundary arguments are connected. It is not a near-wall
turbulence model or external PeleC/physical validation.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, existing fixed tests-enabled configuration.
`test_spray_core` built. Direct `les_dynamic_ghosted` test: **1 passed, exit 0**,
no failures or skips. It verifies:

- Existing exact periodic RHS/face parity and short-halo/invalid-gas rollback.
- Adiabatic thermal closure preserves the uniform-temperature work flux.
- Hot isothermal planes add energy while preserving momentum flux exactly.
- Returned face divergence includes wall heat; doubling the wall temperature
  difference doubles the heat contribution (relative bound `1e-12`).
- Invalid wall temperatures preserve all ghosted RHS/face outputs.

The thermal tests exercise the kernel with controlled donor fields; they do not
claim an integrated reflected-wall AMR solution. No existing tolerance changed.
Provenance preflight, Python syntax, project contract and whitespace passed.
No MPI run, selected build or full configuration matrix was added for this
local kernel change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
build/spray-amr-reactive-debug/tests/test_spray_core les_dynamic_ghosted
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The executable is `build/spray-amr-reactive-debug/tests/test_spray_core`;
the direct test reports `les_dynamic_ghosted: PASS` on stdout and writes no
result artifact.
