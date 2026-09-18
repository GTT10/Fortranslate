# Serial dynamic LES application — 2026-09-18

Implementation `fe7279c37e5672275a62e138a408d333edf8fc28`, based on `c9c9c8c`.
The shared fixed/selected serial spray driver accepts `les_model='dynamic'`,
`les_test_filter` (box default), `les_test_ratio` (two default), and
`les_grid_ratio` (one default). It forwards them through both pure-liquid and
mixture coupled paths. Dynamic settings enter checkpoint identity; nondefault
filter settings require the dynamic model. Unsupported nonperiodic faces and
invalid filters reject before output. Existing default Smagorinsky checkpoint
context is unchanged. Positive `smagorinsky_constant` enables dynamic LES;
its magnitude does not scale the computed dynamic coefficients. Zero disables
SGS. The existing serial coefficient bound of 0 to 0.5 remains.

The periodic-only library limits, frozen-coefficient timestep heuristic and
regularization/unit caveats remain. MPI application support, general boundaries,
AMR dynamic LES and physical qualification remain open.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; unchanged configurations.
Fixed tests-enabled `build/spray-amr-reactive-debug` built `pelef_spray_3d`.
Selected H2O2 tests-disabled `build/spray-amr-reactive-selected-debug` built
`pelef_spray_3d_selected`.

- Existing fixed Smagorinsky boundary/restart CTest passed, 6.60 seconds.
- Fixed dynamic application CTest passed, 2.00 seconds (2.04 seconds total).
- The same dynamic checker directly exercised the selected executable and passed.
- Dynamic fixture: periodic 4-cubed grid, shear, dilute inert drag-coupled parcels,
  Gaussian test ratio four, grid ratio two. Continuous and stop/resume gas and
  parcel CSVs match byte-for-byte; final history matches exactly. Existing mass,
  momentum and energy budget gates remain unchanged.
- Changed/removed model, filter type, test ratio and grid ratio reject restart;
  invalid ratio and nonperiodic faces reject input before history output.
- Provenance, checker/contract Python syntax, project contract and diff whitespace
  checks passed. No skipped tests. No full matrix, MPI dynamic runtime, external
  solver parity or physical validation ran.

Initial CTest passed the existing test but failed the dynamic checker after
successful CSV/history replay because periodic history omits `boundary_mass`.
The checker was extended to use zero exterior exchange for periodic cases and
to assert that periodic history contains no boundary columns. Only the failed
dynamic CTest was rerun, then the selected checker was run directly. No solver
change or tolerance change was needed for that failure.

## Commands and artifacts

Logs: `build/les-dynamic-application-checks/{provenance,build,build-selected,ctest,ctest-dynamic-final,selected-runtime}.log`.
Fixed fixtures are under `build/spray-amr-reactive-debug/tests/spray_dynamic_boundaries/`.
Selected fixtures are under `build/les-dynamic-application-checks/selected-runtime/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(dynamic_boundary_restart|sgs_boundary_restart)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_dynamic_boundary_restart$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tools/check_spray_boundaries.py --sgs .17 --les-model dynamic --executable build/spray-amr-reactive-selected-debug/pelef_spray_3d_selected --work-dir build/les-dynamic-application-checks/selected-runtime
python3 -m py_compile tools/check_spray_boundaries.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
