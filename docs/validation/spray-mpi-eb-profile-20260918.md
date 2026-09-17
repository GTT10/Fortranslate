# MPI EB port mass-flow profiles — 2026-09-18

Implementation `244c6b8a6a03fdf481e348b0c567b3b118a3795a`, base `0425910`.

Port-specific 2--64 point time/rate tables now feed the existing
`integrate_spray_injection_profile` routine. Each selected event explicitly
supplies its integration end time. The integrated area times its event scale
replaces the prescribed packet mass. The template still sets drop size,
composition and temperature; multiplicity supplies the integrated amount.
Zero integrated mass consumes no particle IDs, and restart reconstructs its
expected watermark from effective positive events. All rate points and interval
ends enter the normalized checkpoint context. Single-port/prescribed-mass inputs
remain the default when profile point counts are zero.

This is packet injection at interval start, not continuous mass addition at
every gas step. Intervals may overlap if specified; events are not generated
from the profile automatically. Rates retain the existing 0--1 kg/s bounds.
Integrated mass is limited to 1e6 kg and a 1e6 template-mass factor. Profiles are
namelist arrays; external file loading is not implemented by this change.

## Completed focused checks

GNU Fortran 13.3.0 `/usr/bin/gfortran`, Open MPI 4.1.6 `/usr/bin/mpifort` and
`/usr/bin/mpiexec`, Debug/Ninja. Existing fixed tests-enabled and selected H2/O2
tests-disabled configurations were reused unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_profiles$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-profile-checks/selected-final --selected --scheduled --restart --regrid --ports --profiles
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, test **43.62 s**, total **43.66 s**. Selected runtime script passed. Provenance,
Python syntax, contract and whitespace passed. No failures, skips or tolerance
changes. An initial passing fixture used an explicit zero event multiplier;
the final fixture changes that multiplier to one and uses a zero-width interval,
so the zero-mass result specifically exercises integrated-mass ID handling.
These are revisions of the same test, not extra distinct coverage.

Port one uses constant 0.02 kg/s over 0--40 ns. Port two uses a triangular
0/0.08/0 kg/s profile at 0/20/40 ns. Event start times are [5,5,5,10,20,30] ns,
ends [10,5,10,15,25,35] ns, ports [1,2,2,1,2,1], scales [1,1,1,.5,2,1], and the
run ends at 20 ns. The expected processed total is 1e-9 kg, checked within
relative 1e-14. Ten consecutive IDs and correct port/component state show that
the zero-area event reserves none. Final-time parcels remain at their injection
state. The underlying two-port fixture uses synthetic water/H2 constant-property
tables; it is not physical liquid-hydrogen validation.

Two/four-rank output comparisons retain relative 3e-12, floor 1e-30. A two-rank
root-only checkpoint to four-rank continuation yields exactly equal parsed gas
and particle CSV records to uninterrupted execution after subsequent refinement.
Initial/final and accumulated ledger comparisons retain relative 3e-12. Existing
context, overwrite and truncated-payload rejection checks also run. The fixed
case uses elementary gas initialization; selected gas is pure N2. Chemistry is
implicit and all molecular transport terms are enabled.

Artifacts: `build/spray-mpi-eb-profile-checks/{build.log,build-selected.log,ctest-final.log,selected-final.log}`.
Selected final cases: `selected-final/eb_spray_az1ds5lc/`; CTest LastTest.log records
the fixed directory. Initial successful fixture logs are retained as well.
No full matrix, external solver parity, scaling, convergence, installed-binary
audit or physical validation was run in this increment.
