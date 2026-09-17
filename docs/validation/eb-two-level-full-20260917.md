# Combined two-level EB reaction/transport/hydro — 2026-09-17

Implementation `55fb9463b5058ded10a726abba867f95d30cf320`, based on
`b8ff924` on local `main`.

## Behavior

`advance_reactive_eb_two_level_full_3d` connects reaction half-step, synchronized
transport half-step, synchronized hydro step, transport half-step and reaction
half-step. It reuses the existing masked chemistry, parent donors, aperture flux
restriction, reflux and conservative averaging. The original reaction/hydro API
forwards through the shared split implementation with transport absent.

The operation accumulates physical-boundary flux inventories from hydro and both
transport stages, retains the hydro EB wall contribution and measures chemical
inventory changes across the reaction stages. The minimum species-limiter factor
is the minimum returned by both transport stages. All state arrays and all
inventory/limiter outputs commit together. A later failure discards successful
earlier stages as part of the same transaction.

Disabled transport follows the existing reaction/hydro path and returns limiter
one. Disabled chemistry skips both reaction stages and returns zero chemical
inventory. With both disabled, the existing hydro step is used directly.
The symmetric arrangement does not upgrade the first-order Euler component
steps to an overall second-order method.

Support remains one static serial fine patch, common timesteps, PCM donors,
outflow hydro boundaries and zero physical/EB molecular wall flux. Dynamic or
multiple patches, subcycling, application/checkpoint, MPI, convergence and
physical validation remain open. This is a combined library step, not an overall
project-completion claim.

## Independently completed checks

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; matching existing build
`build/selected-tree-mpi-debug`, MPI on, spray off. Core and selected EB runtime
compiled; execution used fixed elementary chemistry. No MPI or application matrix
was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_two_level_chemistry_3d test_eb_two_level_transport_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_two_level_full_3d test_eb_two_level_chemistry_3d test_eb_two_level_transport_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_two_level_full_3d|unit_eb_two_level_chemistry_3d|unit_eb_two_level_transport_3d)$' --output-on-failure
# After adding the chemistry-disabled/transport-enabled branch assertion:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_two_level_full_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_two_level_full_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

All **3 relevant tests passed**, with no failed or skipped checks:

- Existing reaction/hydro: 2.91 s.
- Existing synchronized transport: 0.49 s.
- Final combined-physics unit: 8.53 s (CTest 1/1, total 8.56 s).
  The earlier three-test run passed in 11.92 s; its 8.48 s combined-unit result
  is superseded by the final rerun and is not counted as an additional test.

The combined unit advances a cut-plane 4-cubed coarse domain and ratio-two
4-cubed fine patch by 1e-8 s with explicit and implicit chemistry. It checks
nonzero reaction and transport responses, positive temperatures, component-wise
boundary/wall/chemical inventory balance, chemical mass/momentum/energy and H/O/N
conservation, inactive refined parents and unchanged covered storage. Disabled
transport matches the reaction/hydro reference exactly; both disabled match the
hydro reference. Chemistry-disabled transport/hydro also satisfies its inventory
balance. Invalid chemistry selection and hydro failure after preceding stages
preserve all four state arrays, three inventories and the limiter output.

Provenance and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-two-level-full-checks/{build,final-build,test-build,ctest,final-ctest}.log`.
External PeleC parity and general EB-AMR qualification are not established here.
