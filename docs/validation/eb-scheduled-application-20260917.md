# Scheduled EB application regridding — 2026-09-17

Implementation: `7c3b9abe4e5012a9d817eee63971d79c33861198`, based on
`b878d21`. The fixed and configure-time selected two-level applications now
connect the conservative remapper to prescribed events and checkpoint restart.

## Implemented behavior

The namelist accepts up to 16 ordered positive event times, each with 1–16
nonoverlapping coarse-index fine-patch boxes. All events are validated before
initialization; events may lie beyond the current final time. The finest
geometry lattice and refinement ratio remain fixed. Timesteps clip to the next
event, including with automatic timestep selection and fine subcycling.

An accepted step landing on an event is followed by conservative remapping,
composite-budget validation, then checkpoint/stop handling. Patch pointers,
coarse coverage and output topology are updated with the new layout. The
original initial inventories and boundary/wall/chemical budgets persist.

Scheduled checkpoints use a separate V1/V2/V3/V4 envelope for fixed/automatic
common stepping and fixed/automatic subcycling. They retain original topology,
all event times/boxes/counts and the completed-event index. Restart verifies
this policy, constructs the expected destination layout, validates all saved
fields and budgets, and checks event index against saved time before publishing
restored solution fields. Inconsistent policy/stage rejects before output.
Static checkpoint layouts remain unchanged. Storage-failure atomic replacement
is not added.

## Focused evidence

GNU Fortran 13.3.0, `/usr/bin/gfortran`, Debug/Ninja, existing
`build/selected-tree-mpi-debug`. Fixed elementary H2/O2 uses explicit chemistry;
selected full H2/O2 uses implicit chemistry with H2/O2/N2 mole fractions
0.29572/0.14785/0.55643.

Both cases use an 8-cubed finest lattice, ratio 2, a planar wall at z=0.31,
reaction/transport/hydro, and two events at 1.5e-8 and 2.5e-8 s. Patch counts
change 2→1→2, merging initial adjacent patches and then moving to separate
boxes. Four accepted steps reach 3e-8 s with two completed regrids. Final
output has 160 active cells, including 64 cells in each fine patch.

| Runtime / timestep mode | Composite budget error |
| --- | --- |
| Fixed / fixed common step | 3.5860547634892602e-16 |
| Fixed / automatic subcycling | 2.2407301912036123e-16 |
| Selected / fixed common step | 4.9316741076790803e-16 |
| Selected / automatic subcycling | 4.6378101363360439e-16 |

For each combination, checkpoints after accepted steps 1, 2 and 3 exercise
restart before the first event, exactly after the first event, and exactly
after the second event. Final solution CSV and cumulative budget CSV are
byte-identical to uninterrupted execution. Changed schedules and corrupted
event indices reject while preserving previous outputs. Existing static
fixed-runtime application regression also passes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_two_level_3d test_selected_full_reactive_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(regression_(selected_)?scheduled_eb_two_level_3d|regression_reactive_eb_two_level_3d)$' --output-on-failure --output-log build/eb-scheduled-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_selected_scheduled_eb_two_level_3d$' --output-on-failure --output-log build/eb-scheduled-checks/selected-ctest.log
python3 -m py_compile tools/check_scheduled_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The first grouped run passed static fixed (7.66 s) and scheduled fixed
(8.22 s), but selected failed because its reused input omitted required
composition. A dedicated selected input was added; its final rerun passed
**1/1**, 18.24 s (total 18.27 s). Thus all three distinct targeted tests have
completed passes; the initial failed run is retained. No skips. Provenance,
Python syntax, project contract and whitespace checks passed.

Build/test logs: `build/eb-scheduled-checks/`. Per-invocation inputs, logs,
checkpoints and CSVs are under
`build/selected-tree-mpi-debug/tests/scheduled_eb_two_level_3d/` and
`selected_scheduled_eb_two_level_3d/`.

This is serial scheduled two-level EB integration. Temperature tagging,
application ratio changes or zero fine patches, general multilevel/MPI EB,
nonplanar schedule qualification, external PeleC parity and physical/convergence
validation remain open. No full matrix or installed-binary audit was rerun.
