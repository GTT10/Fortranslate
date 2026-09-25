# Two-level AMR coefficient input and restart — 2026-09-25

## Source and change

Code `a67a8f7`, parent `104cd9d`. Fixed/selected serial and MPI two-level AMR
front ends enable `difmag`; shared application drivers forward its resolved
value into reacting/transport stepping. Automatic -1 means 0.1 for PPM/pelec
and zero otherwise. Explicit nonnegative input requires PPM.

The existing checkpoint reconstruction text record now contains
`characteristic_ppm difmag-v1 <resolved coefficient>` for PPM. The writer and
both reader paths compare the same canonical representation. Equivalent
automatic/explicit coefficients share identity; changed coefficients reject.
Direct checkpoint configuration validation also rejects nonfinite/negative
coefficients or a positive coefficient outside PPM. The file header and record
count are unchanged. Old PPM records without coefficient identity reject;
PCM/PLM reconstruction records remain unchanged. State payloads are unchanged.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3. Fixed tests-enabled
`build/spray-amr-reactive-debug`; selected H2/O2 tests-disabled
`build/spray-amr-reactive-selected-debug`.

Fixed CTest **1/1 passed**, 47.93 s (47.97 s elapsed). Selected standalone
checker **PASS**. Both serial/MPI executable pairs built. Provenance, Python
syntax, project contract and whitespace checks passed. No failures or skips
in this increment. Compiler warnings remain; selected MPI build dependency
audit passed (ABI 40).

The reacting pressure-outlet fixture includes molecular transport, physical
walls and regridding. Four-rank checkpoint creation followed by two-rank
continuation produces coarse/fine CSVs exactly equal to the serial reference.
Automatic -1 resumed from explicit 0.1 also produces identical CSVs. Coefficient
0.2, changed solver/scheme, wall/pressure changes, chemistry tolerance, integrator,
reaction-rate and mechanism-identity changes reject without final coarse output.
Chemistry-disabled output differs, ensuring the fixture exercises chemistry.
The boundary-ledger comparison retains relative 2e-11 plus absolute 1e-25,
with a shared momentum scale. No tolerance was relaxed.

This is bounded fixed/selected H2/O2 application/restart evidence, not an
installed audit, full matrix, external PeleC oracle or physical validation.
PCM/PLM and EB artificial-viscosity dispatch remain outstanding work.

## Commands and artifacts

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_amr_reactive_3d pelef_mpi_amr_reactive_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_amr_reactive_3d_selected pelef_mpi_amr_reactive_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^mpi_amr_two_level_difmag_restart_3d$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_amr_physical_boundary_restart.py --selected --ppm --difmag --riemann-solver pelec --pressure-outlet --reacting --executable build/spray-amr-reactive-selected-debug/pelef_mpi_amr_reactive_3d_selected --reference-executable build/spray-amr-reactive-selected-debug/pelef_amr_reactive_3d_selected --mpiexec /usr/bin/mpiexec --ranks 4 --restart-ranks 2 --work-dir build/difmag-two-level-input-checks/selected
python3 -m py_compile tools/check_amr_physical_boundary_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/difmag-two-level-input-checks/` contains build, selected build,
CTest, copied `LastTest.log` and selected runtime logs. Fixed fixture path is
recorded in `LastTest.log`. Selected fixture:
`build/difmag-two-level-input-checks/selected/run-8qrzau_p/`.
