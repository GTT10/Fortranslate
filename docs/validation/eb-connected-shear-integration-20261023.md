# Connected-direction production-path shear comparison — 2026-10-23

## Purpose and acceptance

The application lifecycle integration passed nine tests on 2026-10-22.
This follow-up checks whether the production serial recursive solver path
reproduces the frozen isolated candidate's analytic-shear results, using the
new explicit policy argument instead of replacing a library object.
It does not establish local-error convergence or general-geometry accuracy.

`tools/reference_probes/eb_tree_reacting_probe.F90` now accepts
`--audit-shear-connected`. It selects the existing four shear policies and
enables connected directions only for the three separate-momentum policies.
The default policy remains unchanged. The original `--audit-shear` interface
keeps its prior behavior. Both common and recursive native probe calls forward
the flag; this comparison uses the saved recursive case.

The bounded comparison reuses frozen 16-grid, 50 ns cut/regular inputs from
`build/eb-shear-tangent-candidate-20261022/`. It runs the production-linked
probe once per geometry, then compares every CSV row for native diagnostics,
component inventories, all-variable L1/L2/Linf errors and maximum-error
locations. Three enabled policies must exactly match the frozen candidate;
the default must exactly match the original localization results through the
same first step. The regular baseline uses the recovered regular geometry.
No original candidate or 8/16/32 convergence calculation is relaunched.

## Completed results

Status: **PASS**, build and both geometry runs returned zero. Cut took
52.77 s; regular took 102.17 s. No failed or unfinished runs remain.

Independent `audit.py` reread all CSVs, verified source/input/binary and
baseline/output hashes, reconstructed component budget residuals and checked
finite ordered norms plus their matching maximum-error locations.
All **32 CSV tables match exactly at the parsed field level**: 24 enabled-policy
tables match the frozen candidate; eight default tables match the original
results through the first step. These contain 224 inventory rows, 272
17-variable norm rows and 272 maximum-location rows. This establishes equality
of recorded diagnostics, not a comparison of every saved solution cell.

| Geometry | Default budget residual | Momentum / velocity / paired residual |
| --- | ---: | ---: |
| cut | 5.92416e-14 | 5.70231e-14 |
| regular | 6.05880e-14 | 5.60061e-14 |

All correction-accounted normalized budgets pass the existing 1e-12 gate.
The three enabled policies reproduce the frozen candidate's all-variable
L1/L2/Linf errors and maximum-error cell locations. The preserved default
avoids the candidate's previously observed default-policy degradation.
The equality gate and independent audit are saved in `record.json` and
`audit.json`.
Artifacts: `build/eb-connected-shear-integration-20261023/`.
Driver PID at launch: `2251514`. `record.json` stores source/input/executable
hashes, exact commands, return codes, per-file comparison outcomes and budgets.
`run.py` built `test_eb_tree_reacting_reference` in the existing Debug GNU
Fortran/OpenMPI selected-runtime build, then completed the two geometry runs.
`build.log`, `cut.log`, `regular.log` and output CSVs are retained.

The driver saved terminal status and successfully queued its completion
message (return code zero). No calculation was restarted during this audit;
no recovery is required. The equality assertions and independent checks all
passed without changing tolerances.

Provenance preflight, Python syntax, project contract and `git diff --check`
passed during setup and the terminal evidence update (including audit-script
Python syntax). This comparison connects the frozen candidate diagnostics to
the production serial recursive solver at 16 cells and 50 ns. It does not
extend that production-path comparison to 8/32 cells, 200 ns or MPI shear.
The prior study's nonmonotone transverse maxima and weak local convergence
remain open.
See the [application record](eb-connected-application-20261022.md) and
[analytic shear evidence](eb-shear-analytic-20261022.md).
