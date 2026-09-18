# MPI dynamic LES spray application — 2026-09-18

Implementation `678ea2ec6c2e019cfb21ef8baf3bb22048928057`, based on `a074b0d`.
`advance_mpi_mixture_reaction` accepts `les_model='dynamic'` with trailing
`les_test_filter`, `les_test_ratio`, `les_grid_ratio` options (box/two/one by
default). Filter model and ratios join collective control agreement before any
candidate mutation. Dynamic SGS half steps surround hydro inside molecular
transport; coefficients and state are distributed. Periodic faces are required,
even when SGS is disabled. Positive `les_settings(1)` enables dynamic SGS but its
magnitude and fixed Prandtl/Schmidt settings do not determine dynamic coefficients.
Zero retains the disabled convention. The existing candidate gas/particle/film/
ID transaction is retained.

Fixed/selected MPI applications expose the same controls in `mpi_spray` and bind
them to checkpoint identity. Nondefault filter controls require the dynamic model.
Changed/removed settings reject restart; invalid filters and nonperiodic faces
reject before output. The standalone algebraic `advance_mpi_spray_les` API remains
for Smagorinsky/WALE/Vreman; dynamic coupled dispatch uses its dedicated integrator.
Physical boundaries, AMR dynamic LES, large-scale performance and physical
qualification remain open. Existing frozen-coefficient timestep and numerical
regularizer caveats apply.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; unchanged configurations.
Fixed tests-enabled `build/spray-amr-reactive-debug` built `pelef_mpi_spray_3d`.
Selected H2O2 tests-disabled `build/spray-amr-reactive-selected-debug` built
`pelef_mpi_spray_3d_selected`.

- Fixed focused CTest **2/2 passed**, 14.68 seconds, no failures/skips: existing
  WALE restart and new dynamic restart.
- The selected executable passed the same dynamic checker directly.
- Dynamic fixture: periodic 8-by-4-by-4 reacting mixed spray, Gaussian test ratio
  four/grid ratio two, two steps of `1e-8`. Same-rank stop/resume gas/particle
  records match exactly; 2-to-4-rank restart passes existing relative `2e-11`,
  absolute `1e-25` gates. Existing mass/energy, heating, component inventory and
  overwrite-protection checks pass.
- Changed/removed model, filter, test/grid ratios reject checkpoint restore;
  invalid ratio and nonperiodic boundaries reject input. The fixture starts
  with uniform gas; prescribed-gradient MPI coefficient/flux/step evidence
  remains in the preceding records.
- Provenance, checker/contract Python syntax, project contract and whitespace
  checks passed. No full matrix, external solver parity or physical validation ran.

## Commands and artifacts

Logs: `build/les-mpi-dynamic-application-checks/{provenance,build,build-selected,ctest,selected-runtime}.log`.
Fixed fixture artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_dynamic_restart/`.
Selected artifacts: `build/les-mpi-dynamic-application-checks/selected-runtime/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(dynamic|wale)_restart$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_mpi_spray_application.py --executable build/spray-amr-reactive-selected-debug/pelef_mpi_spray_3d_selected --mpiexec /usr/bin/mpiexec --ranks 2 --input cases/mpi_spray_3d/periodic_mixture.nml --work-dir build/les-mpi-dynamic-application-checks/selected-runtime --les-model dynamic
python3 -m py_compile tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
