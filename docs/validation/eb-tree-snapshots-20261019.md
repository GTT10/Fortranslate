# EB-tree intermediate gas output, 2026-10-19

The gas EB-tree application previously wrote only the final composite field.
`output_interval_steps` now enables initial/restored, periodic accepted-step,
and final/checkpoint-stop snapshots. Each snapshot carries the existing gas
CSV and cumulative budget plus a physical-time/step/topology metadata CSV.
Output follows synchronization, rebalancing and regridding. Rejected trial
states are not published. The default remains final-only; this is a gas
output feature, not particle/film output or fixed-PeleC trajectory evidence.
Usage and filenames are in the [case guide](../../cases/reactive_eb_tree_3d/README.md#intermediate-gas-output).

## Completed evidence

The existing `regression_selected_eb_tree_hot_cvode_restart` check passed in
21.85 s, with five physical executions and two expected input refusals.
Twelve snapshots across serial, two-rank, checkpoint-stop and three-rank
continuation runs passed independent gas/species/element/UEINT audits.
Fields at matching steps agree byte for byte. Output intervals 1 and 2
exercise the off-cadence endpoint; restart changes the interval and retains
the saved absolute step/time. Negative intervals and a lexically aliased
checkpoint path refuse before creating output. The reaction-off run retains
default final-only output.

| Serial step | Physical time (ns) | Nodes after regrid | H2 consumed / initial H2 |
| ---: | ---: | ---: | ---: |
| 0 | 0 | 2 | 0 |
| 1 | 100 | 1 | 2.7035764e-5 |
| 2 | 200 | 2 | 9.0332105e-5 |
| 3 | 300 | 2 | 2.2952262e-4 |

All five final gas CSVs and their budget CSVs also match the preceding
no-snapshot binary's saved results byte for byte. This comparison reuses
the completed [CVODE case](reacting-eb-amr-20261019.md#higher-temperature-cvode-extension);
the preceding execution was not repeated. The old binary rejects the new
namelist setting, reproducing the missing capability. Its first diagnostic
log assertion expected the word `namelist`; the actual message was
`Invalid eb_tree input`. The assertion interpretation was corrected from
the saved log without repeating that run.

## Source and reproduction

Dirty `main` base: `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
Changes are limited to the shared EB-tree application output path, the
existing hot-case checker and documentation. No numerical kernel or solver
tolerance changed. Fixed/selected and serial/MPI gas frontends use the
shared source; runtime evidence here is selected full H2/O2, common stepping.

Artifacts: `build/eb-tree-snapshots-20261019/`. `entry.json`, `before/`,
`source/` and `change.patch` isolate this change from pre-existing work.
`build-source-sha256.json` identifies all numerical/application and CMake
sources and was checked again after execution. `result.json` identifies
the binaries, report, prior-output comparison and limitations. `job.json`
contains exact build/test commands and exit codes; logs and JUnit are beside
it. The completed detailed application report is:

```text
build/reacting-hot-diagnosis-20261019/registration/tests/selected_eb_tree_hot_cvode_restart/reacting_restart_0jkms2sx/checks.json
```

The previously configured, unbuilt `registration` directory was built from
the current working tree with GNU Fortran 13.3.0, OpenMPI 4.1.6, Debug,
MPI/spray/SUNDIALS ON, and pinned SUNDIALS 7.2.0 static double-precision
Fortran bindings. The final application-only incremental build took 4.50 s.
Existing compiler warnings remain in the build log; this is not a warning-free
build or a broad configuration audit.

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C HWLOC_COMPONENTS=-gl,-x86 \
  OMPI_MCA_rmaps_base_oversubscribe=1 \
  ctest --test-dir build/reacting-hot-diagnosis-20261019/registration \
  -R '^regression_selected_eb_tree_hot_cvode_restart$' --output-on-failure
```

Mechanism provenance, Python syntax, project contract and whitespace checks
passed (`static-checks.json`). The CTest run had one pass, no failures and
no skips. The prior common/recursive restart matrix and installed-binary
audit were not rerun. This advances the output needed by roadmap C2;
external-solver comparison and physical ignition validation remain separate.
