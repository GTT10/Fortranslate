# SDC forced-reactor backend bridge — 2026-09-27

Implementation `315adb23`, based on `a4af1ecd`, local main, version 0.246.0.

## Change

The shared `forced_reactor_integrator` interface describes one transactional
forced interval. `react_sdc_state_3d` accepts an optional `evaluate_reaction`
procedure and invokes it for each cell. The regular-grid SDC driver forwards
the procedure on every iteration. Omission preserves the native forced
implicit solver. Supplying a backend without a reaction array is rejected.
No SUNDIALS dependency is introduced into the core solver.

The existing bundle metadata permits `explicit` or `implicit`, not `cvode`.
Thus this change adds a backend selection boundary; it does not repair a
previously accepted bundle CVODE setting. Application input and hierarchy
propagation still need to expose this boundary explicitly.

## Independently completed checks

GNU Fortran 13.3.0, Ninja Debug. CVODE uses SUNDIALS 7.2.0 and the selected
runtime in `build/forced-cvode-debug`, configured as in the
[backend record](forced-cvode-20260927.md). The default native path uses the
existing fixed, SUNDIALS-disabled `build/spray-amr-reactive-debug` build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/forced-cvode-debug --target test_sdc_cvode_bridge --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/forced-cvode-debug -R '^unit_sdc_cvode_bridge$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Both focused tests passed: CVODE bridge 0.20 s (CTest total 0.22 s), native
regression CTest total 24.16 s; no failed/skipped tests in the final runs.
Provenance, builds, Python syntax, project contract and whitespace checks passed.
The first test compile used the primitive-conversion arguments in the wrong
order; correcting the fixture call resolved it. Existing real-comparison and
executable-stack warnings remain. Artifacts are in
`build/sdc-cvode-bridge-checks/`, including both preserved LastTest logs.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The new test uses the existing two-species H/H2 fixture on a periodic 4³ grid,
two SDC iterations and a 1e-10 s interval. Constant species sources [0.1,-0.02]
and internal-energy source 1e6 are nonzero. Exactly 128 backend calls occur.
Every cell's final species densities match the independent single-cell forced
CVODE solve exactly; temperature, external mass/energy inventories and reaction
source history meet the explicit test criteria. Chemistry measurably changes
the state. A deliberately destructive failure in the first cell of iteration
2 leaves the entire caller gas state, temperature and reaction history unchanged.
The existing native PPM/SDC regression also passes without a callback.

## Limits and next integration

This verifies library dispatch, not a user-selectable application mode.
AMR, EB, MPI, spray, checkpoint policy, nonuniform reacting fields, detailed
mechanisms, full configuration matrix, installed binaries, full external
solver parity and physical validation were not tested here.
