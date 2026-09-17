# Fixed/selected MPI AMR spray applications — 2026-09-18

Implementation: `0356fd8363496a7afc591ac15238fefb6d28cb1c`, based on
`b3e34d1`. `pelef_mpi_amr_spray_3d` and its selected-mechanism counterpart share
a driver for periodic owner-local gas initialization, mixed-cone initial
injection with collective IDs, coupled particle/gas steps and one optional
prescribed regrid. The driver uses a finest-grid CFL/transport bound even when
recursive gas steps are selected, and retries failed whole steps at half the
interval. Outputs are exposed-cell gas CSV and per-rank particle CSV, plus an
injection ledger and execution summary on standard output.

Selected frontend generation uses the existing mechanism-bundle builder and
its source/SHA contract. The selected executable link edge contains selected
runtimes and excludes `pelef_core`. The case guide documents supported inputs
and the remaining checkpoint, timed-injection, adaptation and boundary work.

## Completed focused verification

GNU Fortran 13.3.0 Debug/Ninja with `/usr/bin/gfortran`, Open MPI
`/usr/bin/mpifort`, `/usr/bin/mpiexec`. Reused the preceding clean fixed
tests-enabled and selected H2/O2 tests-disabled builds; compiler, MPI and backend
configuration did not change. No install audit was run. Configure commands are
in [the reacting-step evidence](spray-amr-reactive-20260918.md).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_application$' --output-on-failure --output-log build/spray-amr-app-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_amr_spray_3d_selected --mpiexec /usr/bin/mpiexec --source /home/gtt13/projects/pelec --work /home/gtt13/projects/pelec/build/spray-amr-app-checks/selected --selected
python3 -m py_compile tests/application/check_mpi_amr_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed CTest: **1/1 passed**, no skips, **2.74 s** test / **2.78 s** total.
Selected application checker: **passed**. Each checker launches three MPI
runs: a two-rank empty-particle reference, then two- and four-rank mixed
injection cases. Thus six final application runs completed successfully.
All runs cross the prescribed regrid time and finish at 20 ns on the root-only
grid. The checker verifies 64 exposed gas cells, four surviving moving parcels,
temperature evolution, evaporation, total gas-plus-liquid mass and water/peroxide
inventories within 3e-12 of their respective scales, and 2/4-rank gas/particle
CSV agreement at relative tolerance 2e-12 (absolute floor 1e-25).

These CLI fixtures disable gas chemistry and molecular transport; those gas
operators have separate coupled-library evidence. This check does not claim
physical ignition, long-time stability, energy validation of the application,
or a new numerical-order result.

## Development failures and corrections

- The first selected build configured outside the mechanism-bundle variable
  scope and failed. Moving frontend generation into the existing bundle builder
  fixed the missing generated symbols.
- The first fixed run rejected fixture hotspot coordinates outside the small
  domain. The fixture now explicitly places the hotspot inside the domain.
- Both checkers exposed a trailing particle-CSV delimiter. The writer now uses
  a termination-aware separator and explicit component-column headers.
- Both then rejected an incorrect test assumption that liquid temperature must
  rise. Evaporation cooled this fixture. The test now checks temperature
  evolution together with independently bounded mass/component conservation.
  No solver tolerance or validation rule was loosened.

Final provenance, Python syntax, project contract and whitespace checks passed.
Compiler warnings remain. No full matrix, installed-binary audit, restart,
external-solver parity, convergence, scaling or physical spray validation was
performed.

Logs: `build/spray-amr-app-checks/` contains build attempts, fixed test failures,
the final `ctest.log`, and `test-selected.log`. Final fixed run artifacts are in
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_application/{empty,two,four}`;
selected artifacts are in `build/spray-amr-app-checks/selected/{empty,two,four}`.
Each run contains input copies, `run.log`, gas CSV and rank-specific particle CSVs.
