# Static LES isotropic application controls — 2026-09-18

Implementation `18f89e8ed8feff4a8297b76eed5ceb67491c88b8`, based on `24dd594`.

Serial/MPI spray coupling and fixed/selected application namelists accept
`les_ci`, default zero, finite in `[0,1]`. Both static SGS half steps receive it.
MPI coupling includes CI in collective control agreement and enables these
steps when either the eddy coefficient or CI is positive; positive CI requires
the hydro path just like the existing SGS control. Nonzero fixed CI with dynamic
LES is rejected rather than ignored. Nonzero CI is serialized into checkpoint
identity; zero preserves the prior context. Changing or removing nonzero CI
therefore rejects restart. Formula, MPI kernel evidence and stability-estimate
limits remain in the [serial](les-isotropic-20260918.md) and
[MPI](les-mpi-isotropic-20260918.md) records.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug, existing configurations unchanged.
Fixed tests-enabled serial/MPI and selected H2O2 tests-disabled serial/MPI
applications built. Selected executables were compiled, not run.
Focused CTest **4/4 passed**, 28.01 seconds, zero failures/skips:

- `spray_mpi_wale_restart`: 5.79 s
- `spray_wale_boundary_restart`: 7.25 s
- `spray_isotropic_restart`: 7.93 s
- `spray_mpi_isotropic_restart`: 7.00 s

The latter two exercise WALE with CI 0.1: exact same-rank continuation,
MPI 2-to-4-rank restart within existing tolerances, conservation budgets,
changed/removed CI checkpoint rejection without history publication, and
negative CI input rejection. Serial checks use wall/outflow cases with chemistry
disabled; MPI uses the periodic mixed-spray fixture. The former two retain the
zero-CI/default path. No tolerances changed. All-three-model algebra and kernel
parity were checked in the preceding records, not rerun here. This increment
does not independently isolate the coupled CI-only (`cs=0`) evolution in a test.

Provenance, Python syntax, project contract and whitespace checks passed.
No builds/tests failed. No full matrix, external solver parity or physical
validation ran.

## Commands and artifacts

Logs: `build/les-isotropic-application-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d pelef_mpi_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(isotropic_restart|mpi_isotropic_restart|wale_boundary_restart|mpi_wale_restart)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py tools/check_spray_boundaries.py tests/spray/test_mpi_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```
