# Recursive particle MOL connection and checkpoint failure — 2026-09-27

The recursive particle wrapper and application now forward `mol_mode`.
Each level uses one old kick/drift and one new kick; subsequent reaction
iterations reuse that level's new source. Existing physical/virtual/ghost
populations, source correction and post-level synchronization remain active.
MOL omits the CTU predictor-source callback. The default CTU path retains it.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6; existing selected heptane Simple-transport
Debug MPI/spray build, tests/SUNDIALS disabled:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2 > build/mpi-mol-recursive-particle-build-final2.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile build/mpi-mol-recursive-particle/run_checks_final.py
git diff --check
```

These commands passed. The build's first attempt failed on overlong Fortran
lines. A subsequent executable rejected a null procedure-pointer actual under
runtime checking. Explicitly omitting the optional argument in the MOL call
resolved that failure; neither earlier attempt counts as passed.

The application driver `build/mpi-mol-recursive-particle/run_checks_final.py`
is **not passing**. Its first launch completes two `1e-7 s` steps with four
evaporating parcels, wall conduction, no chemistry, root shape 4³, and one
ratio-two child covering coarse cells `[1:2,1:4,1:4]`. The nodal 8³ level set
represents the same `x=0.0007 m` plane as the preceding single-level fixture.
The maximum absolute component budget residual is `2.710505431213761e-20`;
component units differ, so this is not a relative accuracy measure.

The second launch fails while writing the first-step gas checkpoint.
`build/mpi-mol-recursive-particle-final2` retains both launches, and
`build/mpi-mol-recursive-particle/driver-final.log` records failure. Planned
rank-changing restart and mode-mismatch launches were not reached.

GDB replays `snapshot.gdb`, `snapshot2.gdb`, and `negative.gdb` localize the
failure to the per-patch composition integral in checkpoint writing. On root
patch cell `(1,3,3)`, species 1 (fuel) has partial density
`-1.2258530725940247e-43 kg/m³`; the strict nonnegative check rejects it.
`negative.log` records the value and stack. Each replay used a fresh directory
to preserve failed checkpoint artifacts. The earlier null-pointer trace is
retained separately in `trace.log`.

No checkpoint tolerance was relaxed and no negative species value was allowed
through the writer. The origin and conservative repair of this cell-average
undershoot remain open; positive reconstructed faces alone do not prove
positive cell averages after a full recursive step. Reacting particle MOL,
film MOL, rank-changing recursive MOL restart, full-duration reference cases
and upstream temporal parity remain unqualified.
