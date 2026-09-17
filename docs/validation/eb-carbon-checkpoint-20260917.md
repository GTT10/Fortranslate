# Generic elemental EB checkpoint and methanol connection — 2026-09-17

Implementation: `8de71b91d34c6d520e94d71f1eda887546b0fe88`, based on
`c84aef7` on local `main`.

The checkpoint API accepts optional `element_names` and `element_atoms` after
`checkpoint_unit`. Supplying both requires complete selected context and
selects schema 3. The file binds the element count, ordered names and atom
matrix, and stores an element inventory of that size. Missing/invalid context,
changed names/counts/atoms, invalid diagnostics and incomplete records reject
before publication. Existing calls retain fixed schema 1 or selected schema 2
with three inventory components. No tolerance or hash was changed in the
existing checkpoint comparison paths.

The two-level selected frontend now consumes generated element tables, uses
generic patch inventories, and checks composite chemical element changes at
output/checkpoint and before publishing a restart. It supports carbon-bearing
species and inert He without interpreting chemical labels. Its new files use
schema-3 embedded records; old schema-2 two-level selected files are rejected,
not migrated. The single-level frontend retains its previous inventory path.

The new chemical check uses `2e-10 * initial_element_amount` plus an absolute
roundoff bound: `64 * epsilon * max(1,accepted_steps) * initial_species_mass *
max(atoms/molecular_weight)` for each element. This handles zero initial
elements without changing the relative tolerance. Diagnostics print element
changes and scales on failure. The first implementation used a `1e-30` scale
floor; a zero-He fixture exposed its unrealistic `2e-40` kmol bound. Measured
He change was 1.5572111749731056e-17 kmol, while C changed
1.3891009462725504e-19 kmol from 1.0705682092400839e-3 kmol. No solver state,
reaction rates or integrator tolerances were adjusted to remove that roundoff.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing `build/selected-tree-mpi-debug`,
configured with `/usr/bin/mpifort` (OpenMPI 4.1.6); these tests are serial.
The optional fixture uses a previously imported bundle, with its YAML SHA
checked against `references/spray_methanol.json`. It is not downloaded by CMake.

```sh
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/selected-tree-mpi-debug \
  -DPELEF_TEST_EB_METHANOL_BUNDLE=/home/gtt13/projects/pelec/build/port-methanol-reference/ffcm1_red.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_reactive_eb_3d_checkpoint pelef_reactive_eb_two_level_3d \
  test_selected_full_reactive_eb_two_level_3d test_selected_methanol_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(unit_reactive_eb_3d_checkpoint|regression_reactive_eb_two_level_3d|regression_selected_eb_two_level_3d|regression_selected_methanol_eb_two_level_3d)$' \
  --output-on-failure
# After correcting the newly added chemical roundoff bound, rerun affected frontends:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(regression_selected_eb_two_level_3d|regression_selected_methanol_eb_two_level_3d)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_selected_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Latest completed passes across those runs:

| Check | Result | Time |
| --- | --- | --- |
| Legacy checkpoint unit | PASS | 0.05 s |
| Fixed application regression | PASS | 7.78 s |
| H2/O2 selected regression | PASS | 11.43 s |
| Methanol selected regression | PASS | 75.12 s |

The final affected-frontends invocation passed **2/2**, zero failures/skips,
86.57 s total. Earlier runs had a 90 s timeout in the methanol implicit case
and a zero-He roundoff rejection in the explicit case. The final methanol
fixture uses explicit integration at the same timestep/tolerances; H2/O2 still
tests implicit integration. Methanol implicit performance remains unqualified.
An initial configuration attempt pointed at the reference manifest instead of
the imported bundle and was corrected before building.

The 21-species methanol fixture exercises two-patch common stepping, single
fine-patch subcycling and two-patch subcycling. Each completes three steps to
3e-8 s without rejected trials, with nonzero CH3OH chemical change, finite
state, positive density/temperature, species closure and independent CSV
inventory integration. Composite normalized residuals are respectively
3.1642207983781379e-16, 2.3361601128111728e-16 and 1.8894125166048556e-16.
The minimum transport limiter is **0** in all three: this is a connection and
conservation test, not transport accuracy or ignition validation.

All three modes restart to byte-identical state/budget CSVs. Changed element
names/atom counts, bundle SHA, composition and stepping mode reject without
overwriting outputs. A separate corruption probe adds 0.01 to the CH3OH
boundary ledger and subtracts 0.01 from its chemical ledger, preserving the
component balance; it rejects specifically with chemical element imbalance
and leaves both CSVs unchanged. Its artifacts are `chemical_imbalance.*` in
the methanol regression directory. A preceding probe that altered the initial
inventory instead was rejected earlier by the existing initial-state check.

Provenance, syntax, contract and whitespace passed. Methanol YAML SHA256:
`ddbba8c97d929530cfb213da45485cca4f979dccf97c1063721bd855d7b3f259`;
H2/O2 YAML SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-carbon-checkpoint-checks/` (first/second/final CTest logs).
Artifacts: `build/selected-tree-mpi-debug/tests/selected_methanol_eb_two_level_3d/`.
No broad matrix, installed audit, external solver parity or physical validation
was run. MPI, dynamic/multilevel EB and single-level carbon inventory wiring
remain open.
