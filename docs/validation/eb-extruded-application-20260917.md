# Extruded 3D EB application and restart — 2026-09-17

## Source and scope

Base `9fc158b802d0b501e7686cdc216399620ce06203`, local `main`.
Connects section-file geometry to fixed/selected EB3D applications and existing
checkpoint geometry comparison. Default planar configuration and record layout
remain. New mode/axis identity is encoded in the prior axis record. Wall-force
momentum directions are diagnosed as changes; invariant tolerances are unchanged.
Mechanisms and generated source are unchanged.

SHA-256:

- `src/driver/simulation_config_reactive_eb_3d_mod.F90`:
  `a4f0e1413b9b79e13a4e93c153f45f0c21e9c1d592a19e533cc5bf378d8527da`
- `src/driver/reactive_eb_3d_application_mod.F90`:
  `d6408e7368df7b13a6fe73633af0d00cc524b1443cd3cacba0830fc3411b63b4`
- `src/io/reactive_eb_3d_checkpoint_mod.F90`:
  `00d01e2aa47a5466b4175fb81cc9f201cc5148acaf4cfada0488b86f671a4df3`
- `tools/check_reactive_eb_extruded_3d.py`:
  `aa18238429dd10e9b67a2b1d627d6078af2b30a161cf26fa7b71f9e9ded0d850`
- `cases/reactive_eb_extruded_3d/section.dat`:
  `8d2c6a8af251535339191bda5026ebb13566aa66a1c888bb5f1025d4ff53203f`
- `cases/reactive_eb_extruded_3d/reference.nml`:
  `411d3bfb35836c100b01b62382b2ec1204f78cc03061107b64ce2a19d8ac8e64`

## Commands and completed checks

GNU Fortran 13.3.0 Debug, existing `build/selected-tree-mpi-debug` configured with
OpenMPI 4.1.6; these EB application and unit checks execute serially.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_3d test_selected_elementary_reactive_eb_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_3d test_selected_elementary_reactive_eb_3d test_reactive_eb_3d_config test_reactive_eb_3d_checkpoint -j 4
python3 -m py_compile tools/check_reactive_eb_extruded_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(regression_reactive_eb_extruded_3d|unit_reactive_eb_3d_(config|checkpoint))$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Builds passed. Final CTest **3/3 passed**, no failures/skips, 5.66 s total:
configuration 0.00 s, checkpoint 0.04 s, extruded regression 5.58 s. Provenance,
Python syntax, project contract and whitespace checks passed. The intermediate
build's guarded-index warning was removed before the final build/test.

## Evidence

The supplied 10x8x3 circular-section fixture runs reaction, molecular transport
and StateRedist hydro for four steps to 2e-8 s. Fixed/selected elementary outputs
match byte-for-byte. Each frontend stops after step one, reloads the checkpoint
with an identical section at a different path and produces the exact continuous
CSV bytes. A changed valid section rejects restart; wrong section dimensions and
section/output path alias reject without replacing completed output or geometry.
The existing planar configuration and checkpoint unit suites also pass.

The fixed continuous log reports maximum invariant error 4.6535150436085246e-12,
maximum elemental error 1.7635305184306054e-12 and transport theta 1. Density and
temperature remain positive. These are small-fixture implementation checks, not
physical validation or a spatial/temporal convergence result.

Logs: `build/eb-extruded-app-checks/{build,rebuild,test-build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/reactive_eb_extruded_3d/{fixed,selected}/`.

## Limits

Runtime coverage uses z extrusion and elementary chemistry; other extrusion axes
have preceding geometry/redistribution unit evidence, not new application runs.
Section shape is piecewise linear and invariant along one coordinate. No general
3D surface builder, multi-cell agglomeration, EB-AMR-MPI integration, force ledger,
full matrix, installed audit, external solver parity or physical validation was
run. Original plane executables cannot read the new extruded identity.
