# Connected-direction application integration — 2026-10-22

## Scope and current verdict

`connected_directions` is wired into the gas EB-tree application and serial/MPI
common and recursive held-MOL paths. It defaults to false and requires
`separate_momentum`. This record covers application integration, not physical
qualification. The focused integration build/test batch **passed all 9 tests**, with zero
failures, disabled tests or skips. Build and CTest both returned zero.
JUnit and the application records were independently reconciled; all recorded
source hashes still match the checkout (`audit.json`).

The planner preserves the existing all-connected rule when omitted/false.
When enabled it retains connected axis pairs, uses zero slope on blocked
axes and falls back to PCM for disconnected diagonals in the retained
subspace. Existing limiter and child offsets are retained.

Application admission retains the separate-momentum requirements: held MOL,
SDC, no redistribution, no LES and no spray. The option is passed through
hydro/transport halo construction, common/recursive stages, history identity
and replay. MPI selection uses bit 8 and checks collective agreement.
Checkpoints prepend `PELEF_EB_TREE_CONNECTED_DIRECTIONS_V1` only when enabled;
readers reject policy changes in either direction. False does not add a token.

## Reproducibility

Artifacts: `build/eb-connected-application-20261022/`.
`record.json` captures the dirty-main HEAD, per-source SHA256, exact command
arrays, PID, return codes and terminal status. `before/` preserves the
production-file versions before this integration. `verify.py` is the batch
driver; `build-integration.log`, `ctest.log` and `ctest.xml` hold its evidence.

Build directory:
`build/eb-eos-application-20261021/isolated-openmpi/registration`.
GNU Fortran 13.3.0, Debug, configured MPI compiler `/usr/bin/mpifort.openmpi`
and launcher `/usr/bin/mpiexec.openmpi`. The shell-default Intel MPI launcher
is not the launcher selected for these tests.

The initial fixed serial/MPI application build completed successfully
(`build.log`). The integration batch successfully built those applications, selected
full-mechanism serial/MPI applications, MPI halo and serial transfer unit
executables, then ran nine focused CTests sequentially (415.29 s):

- Fixed common/recursive connected-direction application tests.
- Selected hot-CVODE common/recursive connected-direction restart tests,
  combined with paired pressure and converged pinned EOS.
- Existing fixed held and paired-pressure common application controls.
- Two/three-rank halo tests, including collective rejection without output mutation.
- Serial reactive-transfer unit test, including connected stencil admission.

The application checkers exercise serial/MPI agreement, changed-rank restart,
policy mismatch refusal and invalid admission. Both fixed connected-direction checks passed (28 subprocess cases each),
including exact serial/MPI/restart gas, temperature and UEINT comparison.
Maximum reported inventory residuals were 1.097e-15 (common) and 1.132e-15
(recursive). Both selected hot-CVODE checks passed, including 2-to-3-rank
restart, intermediate snapshots, reaction-off control, bidirectional policy
mismatch and mixed-rank refusal. Their maximum recorded budget residuals
were 9.8192e-16 (common, including reaction-off) and 8.1332e-16 (recursive).
Hydrogen reacted fractions were 2.29522552e-4 and 2.29522631e-4; these are
short-time chemistry integration checks, not ignition validation.
Detailed records are under the build directory `tests/eb_tree_connected_application_*`
and `tests/selected_eb_tree_connected_restart_*/reacting_restart_*/checks.json`.
Dedicated connected-direction retry/history corruption testing is not added
by this batch.

## Completed static checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: PASS.
- Python compilation of both modified application checkers and `verify.py`: PASS.
- `python3 tools/check_project_contract.py --source .`: PASS.
- `git diff --check`: PASS after the terminal evidence update.

## Notification and recovery

Driver PID at launch: `2208793`; the batch completed at 23:24:12 JST.
The driver saved terminal status and successfully queued its completion
message (return code 0), confirmed in `completion-notification.json`. The first shell-background launch produced
no driver process or record; the detached Python launch above completed the batch. No duplicate run was
started during this evidence audit.

No failed or unfinished verification remains in this batch. No recovery run
is needed. Further physical convergence work remains a separate gate.

## Scientific boundary

The prior frozen-candidate 8/16/32, 200 ns shear study improved integrated
errors with separate momentum, but transverse local-error convergence remains
open. That isolated executable is distinct from this production integration.
Sphere/corner donor tests establish topology admission, not curved-geometry
flow accuracy or arbitrary cut-cell volume-weighted conservation. Neither
application lifecycle checks nor these candidate results establish PeleC
parity or complete physical qualification. See the
[analytic-shear record](eb-shear-analytic-20261022.md).
