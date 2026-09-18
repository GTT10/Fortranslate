# WALE/Vreman LES grid flux and integration — 2026-09-18

Implementation `dd84758`, base `52bf72d`.

`compute_les_rhs_3d`, `compute_les_ghosted_rhs_3d` and
`advance_les_ssprk2_3d` accept optional trailing `model`, with exact values
`smagorinsky`, `wale` or `vreman`. Omission selects the previous Smagorinsky
path. The existing coefficient argument `cs` supplies Cw for WALE, Cs for the
other two models. Unknown names reject before coefficient/dt-zero shortcuts.

Face velocity gradients select the corresponding deviatoric stress kernel;
its viscosity also drives existing heat/species diffusion and the maximum
momentum/thermal/species diffusion bound. Both RK2 stages use the selected
model. Existing boundary flux accounting and atomic state/temperature publication
are retained. Model selection is available to supplied-halo callers, but MPI
wrappers, spray coupling, CLI inputs and restart identity are not connected yet.
This is not physical validation or full pinned LES solver equivalence.

## Focused checks

GNU Fortran 13.3.0, reused Debug/Ninja fixed tests-enabled and selected H2O2
tests-disabled configurations. No MPI runtime test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_(wale|vreman|grid|zero|reject|boundary)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

New WALE/Vreman tests use crossed sinusoidal velocities on 4^3 cells with
periodic or outflow boundaries. They require positive model diffusivity, all
component flux/boundary balance, integrated boundary budgets, periodic kinetic
energy dissipation, large-step state/temperature rollback and invalid-model
rejection even with zero coefficient/dt. Existing Smagorinsky grid, boundary,
zero-coefficient and rollback cases are retained. Supplied-halo model forwarding
is compiled but not independently runtime-tested here.

CTest **6/6 passed**, zero skips, 0.28 s total. Fixed/selected builds, provenance,
Python syntax, project contract and diff checks passed. No failed or interrupted
runs. Logs: `build/les-alternative-grid-checks/{provenance,build,build-final,build-selected,ctest}.log`.
