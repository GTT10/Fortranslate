# Larger selected mechanism capacity — 2026-09-15

Source: parent `516af42` plus this commit's common capacity, importer/generator,
selected configuration, reactive transport, MPI comparison storage, tests and
documentation. Fortran runtime capacity is now 128 species rather than 32.
The Python generators and CMake target guards use the same limit. This raises
a storage/input bound; it does not qualify every mechanism up to that size.

## Evidence

- The generator suite passed 12 tests in 21.844 seconds, including rejection
  beyond the revised capacity. No reaction/thermo acceptance forms changed.
- GNU Fortran 13.3 Release with SUNDIALS 7.2.0 in the new `build/port-gri53`
  built `pelef0d_selected` and `pelef_spray_3d_selected` for GRI-Mech, 53 species
  and 325 reactions, including its original third-body and falloff reactions.
- Cantera 3.2.0 constant-volume methane/air at 1500 K and 101325 Pa was compared
  with selected CVODE through 0.1 ms. Initial mole amounts CH4:O2:N2 were
  1:2:7.52. Maximum temperature error was 3.5470293369e-11 K; maximum species
  mass-fraction error was 6.7533863828e-15. Methane changed by more than 1e-7
  in mass fraction, preventing an unchanged-state pass.
- A 4x4x4 53-species case enables CVODE chemistry, molecular transport and
  liquid-water drag/heating/evaporation. It runs to 0.1 microseconds, stops
  after one step, and restarts. Final gas/parcel CSVs and history match exactly;
  the energy residual is within 1e-9 of total energy. This is a coupling and
  restart exercise, not a physical methane spray validation.
- `test_reactive_transport_3d` was rebuilt in GNU Debug and run directly;
  it passed. Its invalid species-count probe now uses capacity+1. The existing
  7/10-species transport checks remain unchanged. An initial CTest filter
  matched no tests; it is not counted as a pass.
- Source provenance, project contract, Python syntax and whitespace pass.

Source: Cantera 3.2.0 packaged `data/gri30.yaml`, SHA256
`06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345`.
Source/bundle hashes are recorded in `references/gri53.json`.
Artifacts and machine-readable errors:
`build/port-gri53/check/gri53-0ggaboao/result.json`.

## Commands

Copy the packaged `gri30.yaml` from the Cantera environment into
`build/port-gri53-reference/gri30.yaml`, then:

```sh
build/cantera-venv/bin/python tools/ingest_cantera_mechanism.py --input build/port-gri53-reference/gri30.yaml --output build/port-gri53-reference/gri30.json --phase gri30 --symbol-prefix gri53 --module-name gri53_bundle_mod --chemistry-integrator implicit
cmake -S . -B build/port-gri53 -G Ninja -DCMAKE_BUILD_TYPE=Release -DPELEF_ENABLE_TESTS=OFF -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_SUNDIALS=ON -DCMAKE_PREFIX_PATH=/home/gtt13/projects/pelec/build/sundials-7.2-fortran/install -DPELEF_MECHANISM_BUNDLE=/home/gtt13/projects/pelec/build/port-gri53-reference/gri30.json -DPELEF_MECHANISM_SOURCE=/home/gtt13/projects/pelec/build/port-gri53-reference/gri30.yaml
cmake --build build/port-gri53 --target pelef0d_selected pelef_spray_3d_selected --parallel 8
build/cantera-venv/bin/python tools/check_large_mechanism.py --build-dir build/port-gri53 --source build/port-gri53-reference/gri30.yaml --work-dir build/port-gri53/check
python3 -m unittest discover -s tests/python -p test_generate_mechanism_bundle.py
cmake --build build/port-spray-debug --target test_reactive_transport_3d --parallel 8
build/port-spray-debug/tests/test_reactive_transport_3d
python3 tools/check_mechanism_provenance.py --bundle build/port-gri53-reference/gri30.json --source build/port-gri53-reference/gri30.yaml
python3 tools/check_project_contract.py --source .
python3 -m py_compile tools/check_large_mechanism.py tools/generate_elementary_mechanism.py tools/ingest_cantera_mechanism.py tools/compare_mpi_amr_reactive_3d.py
git diff --check
```

Initial build caught an overlong declaration after enlarging the array name;
it was wrapped. The initial checker omitted a final newline in its namelist
and Fortran rejected parsing; the writer now adds the newline. Subsequent
runs passed without changing numerical tolerances.

No new >32-species MPI/AMR/EB runtime qualification, full compiler matrix,
128-species production benchmark, sparse chemistry solver or diesel mechanism
validation was run. Dense chemistry and diffusion cost still grow with size.
