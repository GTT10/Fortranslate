# Forced chemistry adaptive integration — 2026-09-26

Implementation `80163dad`, parent `7d8e9856`, local main, version 0.246.0.

`advance_forced_reactor_implicit_adaptive` accepts one substep in species-density
variables with fixed external species and internal-energy-density rates.
At each backward-Euler endpoint, density follows its external mass inventory
and specific energy follows rhoe(0) + elapsed * energy_source. The existing
reduced chemical Jacobian and Newton line search solve the endpoint composition.
An optional residual reference separates the externally forced balance from
the physical Newton starting guess; a negative externally predicted species
does not by itself invalidate a reacting trial.

Full-step/two-half-step comparison controls the substep. Richardson
extrapolation is accepted when the forced RHS accepts its state; otherwise
the physical two-half-step result is retained. Failed attempts shrink the
step, and exhaustion preserves caller state and temperature. Absolute
tolerance is in species-density units. The empty-mechanism limit adds only
external sources. Existing constant-volume callers retain their interface
and default residual definition.

## Completed checks

- GNU Fortran 13.3.0 Debug, fixed tests-enabled directory
  `build/spray-amr-reactive-debug`: target `test_implicit_h2o2_reactor` built.
- `unit_implicit_h2o2_reactor`: **1/1 passed**, 0.02 s test / 0.07 s total.
  Added checks cover inert source integration, invalid-step rollback, the
  zero-source limit against the existing adaptive implicit integrator, and
  mass/energy inventories with simultaneous chemistry and external sources.
  Existing RHS and constant-volume checks also remain active.
- Provenance, Python syntax, contract and whitespace passed. No failed or
  skipped tests in this change. Build warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_implicit_h2o2_reactor --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_implicit_h2o2_reactor$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/forced-reactor-integrator-checks/` contains build, provenance,
CTest and LastTest logs. This is a library adaptive substep, not yet SDC
application chemistry. Full-interval subcycling, SDC coupling and external
solver comparison remain open; no physical-validation claim is made.
