# Regular MPI spray difmag controls — 2026-09-25

## Source and implementation

Code commit `dc689f1`, parent `ddd7fe5`. GNU Fortran 13.3.0, OpenMPI,
Ninja Debug and Python 3.12.3. Fixed tests-enabled build:
`build/spray-amr-reactive-debug`; selected H2/O2 tests-disabled build:
`build/spray-amr-reactive-selected-debug`.

`advance_mpi_spray_hydro` and `advance_mpi_mixture_reaction` accept trailing
optional `difmag`. Both resolve omission to 0.1 for PPM/pelec and zero otherwise,
validate finite nonnegative values, and include the coefficient in collective
agreement before state publication. Positive coefficients require PPM.
The coupled operator forwards it into hydro and the final CTU correction.

Fixed/selected MPI spray namelist input accepts automatic `-1` or an explicit
nonnegative PPM coefficient. Checkpoint context contains the resolved PPM value.
Older PPM context is incompatible; default PLM context is unchanged. No binary
payload/schema change was needed.

## Completed checks

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_spray_3d test_spray_mpi_ppm --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_ppm_4|spray_mpi_difmag_restart)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_difmag_restart$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/spray/test_mpi_spray_application.py --executable build/spray-amr-reactive-selected-debug/pelef_mpi_spray_3d_selected --mpiexec /usr/bin/mpiexec --ranks 2 --ppm --difmag .037 --input cases/mpi_spray_3d/periodic_mixture.nml --work-dir build/difmag-mpi-spray-checks/selected
python3 -m py_compile tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

- `spray_mpi_ppm_4`: PASS, 6.60 s. Four two-cell x slabs match serial gas
  state/temperature exactly with coefficients 0, 0.1 and 0.037 across periodic,
  physical wall/outflow and six Mach-inlet orientations. Surface budgets retain
  tolerance 2e-13. Rank-dependent coefficients and one-rank invalid negative
  coefficients reject collectively with gas and temperature unchanged.
- `spray_mpi_difmag_restart`: final PASS, 10.53 s (10.57 s elapsed).
- Selected standalone application checker: final PASS.
- Provenance, syntax, contract and whitespace checks: PASS.

Application checks use coefficient 0.037 in PPM/pelec, preserve exact same-rank
restart snapshots, and retain 2e-11 relative / 1e-25 absolute tolerances for
2-to-4-rank snapshots. Gas species closure and mass/energy budgets retain their
existing bounds. Changed coefficient, -2, NaN and explicit coefficient with PLM
are rejected. Automatic and explicit 0.1 resumed on four ranks have exactly
equal gas/particle snapshots and history rows.

## Failures and correction

Initial builds failed on new lines exceeding Fortran's 132-column limit in
the library, unit test and namelist declaration. Continuations fixed them;
compiler warnings remain. Failed logs are preserved as `first-build.log`,
`second-build.log`, `third-build.log`, `first-selected-build.log` and
`second-selected-build.log`.

The first application checks failed their fixed two-step assertion. The
PPM/pelec coupled case successfully retries the first step at 5e-9, then reaches
1.5e-8 and 2e-8, producing three accepted steps. For this coefficient fixture the
checker now checkpoints at the uninterrupted run's first accepted time, which
preserves the step sequence for exact restart comparison. It checks history
length against the actual final step; older fixtures retain their two-step
requirement. No solver acceptance threshold or numerical comparison tolerance
was changed. The internal reason for the first rejected coupled attempt was
not independently diagnosed here. Initial fixed CTest had 1/2 pass (hydro) and
1/2 fail (application assertion); the corrected application alone was rerun.
Initial selected check failed the same assertion. These are not reported as
an initially clean 2/2 run.

## Artifacts and limits

All logs: `build/difmag-mpi-spray-checks/`. `first-test.log` and
`first-LastTest.log` retain the initial fixed run; `first-selected-test.log`
retains the initial selected assertion. `LastTest.log` points to final fixed
fixture inputs, per-rank CSVs, checkpoints and rejection logs. Final selected
fixture: `build/difmag-mpi-spray-checks/selected/run-ngm1gido/`.

This increment covers regular MPI gas/particle execution. AMR coefficient
propagation/input/restart remains open. No full matrix, installed audit, complete
PeleC oracle, detailed-fuel run or physical spray validation was performed.
