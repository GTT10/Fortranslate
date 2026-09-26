# Coupled serial EB tree SDC — 2026-09-27

## Source and behavior

Implementation `3082ab79`, based on `4ed93b10` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`advance_reactive_eb_tree_sdc_3d` couples source-aware CTU hydro, molecular
transport, external force/work and forced chemistry at a common timestep.
The driver starts from a restricted private tree, evaluates old sources,
iterates with averaged old/new transport and external sources plus the
persistent reaction estimate, and solves forced chemistry from the old state.
Provisional stages retain refined-parent states. Final integrated hydro and
averaged transport fluxes enter the existing peer and coarse/fine
synchronization, including final restriction.

Boundary exchange is calculated from synchronized root face fluxes. Wall
exchange is integrated from the old-state wall source; external exchange uses
the volume integral of time-centered forcing. Chemical exchange integrates
the difference between reacted state and the last nonreacting candidate.
These are separate source calculations, not a closure residual assigned to a
ledger. Tree state, temperature, optional reaction estimates and all diagnostics
publish only after every operation succeeds.

## Independently completed checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- `test_eb_tree_full_3d` built successfully.
- CTest `unit_eb_tree_full_3d`: 1/1 passed, 52.85 s test / 52.90 s total.
- New coupled check: three levels, five patches, z-plane cut at 0.31, warmed
  1000 K fixture, dt 1e-10 s, two SDC iterations, chemistry and all molecular
  transport flags enabled, force [20000, -30000, 40000] per fluid volume.
  Maximum normalized component ledger residual: 3.7650552248896492e-15,
  below the specified 2e-11 bound. Chemistry and forcing ledgers are nonzero.
- Invalid artificial viscosity rejects after successful old transport
  evaluation. Gas state, temperature, reaction estimates, all four exchange
  ledgers and limiter remain exactly unchanged.
- Existing split and reaction-stage tests remain in the same executable.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_full_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-eb-tree-loop-checks/`, including `test-details.log` with the
reported residual. No application trajectory/restart, recursive SDC step,
selected build, MPI run, full matrix, installed-binary audit or physical
validation was performed. This loop uses the current tree's conservative
unredistributed EB rates; single-grid small-cell redistribution options are
not exposed by this tree interface. The bounded conservation check does not
establish temporal convergence or full pinned PeleC parity.
