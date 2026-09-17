# Nodal 3D EB application and restart — 2026-09-17

## Source and implementation

Base `71c58ad7612afc4c21416f7bfdda8508ae32e95a`, local `main`.
Connects the 3D nodal geometry builder to fixed/selected application input,
coupled stepping and geometry-bound restart. Mode `level_set` reads dimensions
and finite nodal values from `level_set_file`, constructs private geometry and
uses the existing metric comparison at restart. The geometry selector stores
`level_set`; old plane/extruded selectors retain their representation. Numerical
kernels, mechanisms, checkpoint tolerances and invariant gates are unchanged.

SHA-256:

- `src/driver/simulation_config_reactive_eb_3d_mod.F90`:
  `0e6fa6ec0f95ab3d5eb9fea988a0e79629fdfc3a1776af6cf18f5bb61e82e96f`
- `src/driver/reactive_eb_3d_application_mod.F90`:
  `24d066a49b6b0fd6a1a23cb354540d4a1fd7799505807587a02682becc926784`
- `src/io/reactive_eb_3d_checkpoint_mod.F90`:
  `0e3e607cc0bd851c6b8dd34e56a5305d7e16dca732eba7e404355fad5f07afff`
- `tools/check_reactive_eb_extruded_3d.py`:
  `19216e8b35245beb81b9d0dbf1bfb263c7a3b3bf7172a2a38c48cac43bc7ccc9`
- `cases/reactive_eb_level_set_3d/reference.nml`:
  `6dd8bc8917a2a643702080c2a1cda7e731f3d1a948b94c75d8d8e648bc4ca3c7`
- `cases/reactive_eb_level_set_3d/level_set.dat`:
  `d614b017359f3e6bb4433c4a87f791a4092a226e868369fc8a8fed2199316be6`

## Commands and completed results

GNU Fortran 13.3.0 Debug in existing `build/selected-tree-mpi-debug` (OpenMPI 4.1.6
configuration). These focused EB unit/application runs execute serially.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_3d test_selected_elementary_reactive_eb_3d test_reactive_eb_3d_config test_reactive_eb_3d_checkpoint -j 4
python3 -m py_compile tools/check_reactive_eb_extruded_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(regression_reactive_eb_level_set_3d|unit_reactive_eb_3d_(config|checkpoint))$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Builds passed. CTest **3/3 passed**, no failures/skips, 12.03 s total:
configuration 0.00 s, checkpoint 0.04 s, nodal application regression 11.96 s.
Provenance, Python syntax, project contract and whitespace checks passed.
Existing compiler warning categories remain.

## Evidence

The supplied 10x8x6 spherical-solid nodal fixture reconstructs 96 cut cells and
16 covered cells. Reaction, viscosity, conduction, species diffusion and
barodiffusion remain enabled. Fixed and selected elementary frontends each
complete four steps to 2e-8 s, save after step one, and resume from an identical
nodal file at a different path. Final output is byte-identical to continuous
execution and between fixed/selected mechanisms. A valid changed sphere rejects
restart; wrong file dimensions and input/output aliasing reject without replacing
completed outputs or geometry input.

The fixed continuous log reports maximum invariant error 3.1132958403061922e-14
and maximum elemental error 2.8337152293529549e-14. Existing plane configuration
and checkpoint units also pass. The previous extruded application branch of the
shared checker was not independently rerun in this increment.

Logs: `build/eb-level-set-app-checks/{build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/reactive_eb_level_set_3d/{fixed,selected}/`.

## Limits

The sphere is a faceted nodal reconstruction. Unresolved/disconnected fluid cells,
degenerate geometry and small-cell clusters without an eligible single receiver
remain unsupported. No detailed-mechanism CFD matrix, convergence, larger-neighborhood
stabilization, EB-AMR-MPI integration, wall-force ledger, full matrix, installed
audit, external solver parity or physical validation was run.
