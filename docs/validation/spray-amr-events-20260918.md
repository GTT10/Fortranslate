# Scheduled MPI AMR spray injection — 2026-09-18

Implementation: `c16c82b08d1b879e5eced431ff7da6a5917b8905`, based on
`cc71fbf`. The shared fixed/selected MPI AMR spray driver now accepts up to
1024 explicit injection times with per-event parcel counts and total masses.
Times are finite, nonnegative and strictly increasing. Count/mass pairs must
both be zero or both positive. With no explicit schedule, the original initial
injection remains the default. Events share the configured cone and mixture.

The timestep is clipped to the next event. Regridding occurs before injection
when their times coincide. Events at final time apply before output; future
events remain unapplied. Injection and finest-owner migration use private
particle/component arrays and a private ID watermark; only their success and
finite cumulative inventory permit publication. The applied-event count is
printed for later restart integration.

## Completed focused verification

Reused GNU Fortran 13.3.0 Debug/Ninja fixed and selected H2/O2 build directories,
with explicit `/usr/bin/gfortran` and Open MPI `/usr/bin/mpifort`,
`/usr/bin/mpiexec`. Compiler/MPI/backend configuration did not change. Fixed
tests were enabled and selected tests disabled. No installed prefix was used.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_application$' --output-on-failure --output-log build/spray-amr-events-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_amr_spray_3d_selected --mpiexec /usr/bin/mpiexec --source /home/gtt13/projects/pelec --work /home/gtt13/projects/pelec/build/spray-amr-events-checks/selected --selected
python3 -m py_compile tests/application/check_mpi_amr_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed CTest **1/1 passed**, no skips, **4.52 s** test / **4.55 s** total.
The selected application checker also passed. Each checker runs an empty
two-rank reference, two/four-rank default initial injection, and two/four-rank
scheduled injection: **ten completed MPI runs** across both executables.
Scheduled events inject four parcels at zero and another four at 10 ns,
coincident with coarsening; final time is 20 ns. Checks confirm two applied
events, eight unique sequential IDs, moving/evaporating parcels, cumulative
gas-plus-liquid mass and water/peroxide inventories within 3e-12 of their
scales, and gas/particle 2/4-rank agreement at relative 2e-12 with absolute
floor 1e-25. Default initial-injection behavior remains covered.

All build/test attempts in this increment passed. Provenance, Python syntax,
project contract and whitespace checks passed. Existing compiler warnings
remain. The fixture disables gas chemistry/transport. No restart, final-time
event-specific test, invalid-schedule CLI test, multiport injection, installed
audit, convergence, scaling or physical spray validation was performed.

Logs: `build/spray-amr-events-checks/{build,build-selected,ctest,test-selected}.log`.
Fixed artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_amr_application/`.
Selected artifacts: `build/spray-amr-events-checks/selected/`.
Each contains `empty`, `two`, `four`, `events_two` and `events_four` runs.
