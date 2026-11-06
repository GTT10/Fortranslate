# Reacting EB-AMR application lifecycle, 2026-10-19

The selected full-H2/O2 application now has a non-inert H2/O2/N2 case through
common and recursive MOL intervals, scheduled child removal/recreation and
changed-rank restart. A focused higher-temperature CVODE extension also
completes regridding and changed-rank restart. This advances completion gate
C1 within short reacting intervals; it does not close stiff ignition, nonzero-material
partial recovery, fixed-PeleC trajectories or production mesh qualification.

## Completed application evidence

`reacting_restart.nml` uses mole proportions 2:1:3.76, initial pressure
150000 Pa, regular/cut density 0.31/0.35 kg/m3, a 4-cubed root, one partial
child, retained internal energy, two SDC iterations and molecular transport.
Three 1e-7 s steps remove the child at step 1 and recreate it at step 2.

| Check | Result | Elapsed seconds |
| --- | --- | ---: |
| `regression_selected_eb_tree_reacting_restart_common` | Passed | 76.3913 |
| `regression_selected_eb_tree_reacting_restart_recursive` | Passed | 104.436 |
| Existing `regression_selected_eb_tree_held_application` pure-N2 control | Passed in the preceding attempt | 54.95 |

Each new application check includes a reaction-off control, serial/2/3-rank
execution, 2-to-3-rank and serial/MPI continuation, exact restore without an
additional step and rejection of changed chemistry. All final gas, temperature,
heat-release and UEINT CSVs agree exactly with uninterrupted execution.
Independent cell integration closes the signed source/boundary budgets;
species balances additionally use their own inventories without a unit-sized
floor. H/O/N chemical inventories and NASA7 mixture internal energy are
checked independently of the driver's printed residual.

Hydrogen consumption is small but resolved: fractions `4.7602650e-9` and
`4.7602673e-9` for common and recursive stepping. The corresponding maximum
serial gas-budget residuals are `7.1445e-16` and `5.4245e-16` using the checker's
documented scale. These values do not establish a strongly reacting regime.

## Source and reproduction

Dirty-main base: `0a0e00ebe9138145fbb0e0e745665d0d768df009`. The complete
frozen source, source manifest, binary hashes, commands, logs and JUnit are
under `build/completion-reacting-amr-20261019`. The numerical snapshot was
checked against its manifest after the completed tests. GNU Fortran 13.3.0,
OpenMPI 4.1.6, CMake 3.28.3, Debug, MPI/spray ON and SUNDIALS OFF were used.

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C HWLOC_COMPONENTS=-gl,-x86 \
  OMPI_MCA_rmaps_base_oversubscribe=1 \
  ctest --test-dir build/completion-reacting-amr-20261019/build \
  -R '^regression_selected_eb_tree_reacting_restart_(common|recursive)$' \
  --output-on-failure --output-junit application-tests.xml
```

The accepted run reports are named in `result.json`; the full fixture has
11 application invocations per clock mode, including expected rejection.
Fourteen provenance regressions, bundle provenance, Python syntax, the
project contract and whitespace checks passed. No broad configuration matrix
or installed-binary audit was run.

The original 30 ns pilot failed the minimum reaction signal. Extending the
case to 300 ns retained the same acceptance threshold. The first full attempt
failed a checker assumption: cross serial/MPI restart retains the saved
initial inventory, while a fresh destination run can group its initial
reduction differently. The corrected check requires the original prefix
exactly and independently audits each species balance. Runtime code and
physical tolerances were not changed for these corrections. Both failed
attempts remain separate in the artifact directory.

## Higher-temperature CVODE extension

`reacting_stiff.nml` raises the initial pressure to 250000 Pa and selects
the existing CVODE backend, retaining the mixture, densities, transport,
mesh lifecycle and 300 ns endpoint above. A bounded native-backend pilot
did not complete a coarse step within its 120 s cost limit. A separate
first-cell debugger observation completed successfully but used 391 internal
steps for a 100 ns interval at 1796.49 K. This identifies substantial local
integration work, not a hang or proof of a numerical defect. No reactor
implementation or reaction tolerance was changed.

The CVODE check passed with five application invocations:

| Invocation | Elapsed seconds |
| --- | ---: |
| Serial, full interval | 6.485 |
| Reaction-off control | 2.772 |
| Two ranks, full interval | 4.881 |
| Two ranks, checkpoint stop | 2.172 |
| Three ranks, checkpoint continuation | 3.226 |

The serial hydrogen consumption fraction was `2.2952261841628982e-4`
(0.02295%), exceeding the unchanged hot-case gate of `1e-5`. Reaction-off
consumption was zero. Serial, two-rank and resumed three-rank final CSVs
were byte-identical. Independent species, element and EOS checks passed;
the maximum scaled serial gas-budget residual was `2.2687e-16`. These
measurements establish neither native/CVODE trajectory parity nor ignition
or production-scale performance. The native pilot is censored, so a complete
runtime speedup ratio cannot be reported.

Artifacts are under `build/reacting-hot-diagnosis-20261019`: `cell.log`
contains the native diagnostic; `configure.log`, `build.log` and
`build-result.json` identify the clean CVODE build and binaries. It uses the
same immutable numerical source snapshot described above, GNU Fortran
13.3.0/OpenMPI 4.1.6 Debug, MPI/spray ON, SUNDIALS ON and the existing pinned
SUNDIALS 7.2.0 static double-precision Fortran installation at
`build/sundials-7.2-fortran/install`. Numerical/application sources were
verified unchanged from the snapshot. The new checker ran from this working
tree, with its hash and the case hash recorded in `build-result.json`.

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C HWLOC_COMPONENTS=-gl,-x86 \
  python3 tests/application/check_eb_tree_reacting_restart.py \
  --serial build/reacting-hot-diagnosis-20261019/cvode-build/tests/test_selected_full_reactive_eb_tree_3d \
  --mpi build/reacting-hot-diagnosis-20261019/cvode-build/tests/test_selected_full_mpi_reactive_eb_tree_3d \
  --mpiexec /usr/bin/mpiexec \
  --work build/reacting-hot-diagnosis-20261019/hot-cvode --hot-cvode
```

The completed report is
`hot-cvode/reacting_restart_06l2u2sa/checks.json`, with inputs, commands,
binary hashes, native logs and CSVs beside it. The new CTest registration
`regression_selected_eb_tree_hot_cvode_restart` was verified by configure
and `ctest -N`; the execution above was a direct checker invocation, not
a CTest execution. Provenance preflight, Python syntax, project contract
and whitespace checks passed. The previously completed broader restart
matrix was not rerun; hot recursive subcycling was not exercised.
