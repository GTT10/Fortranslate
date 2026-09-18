# Characteristic pressure ghost-state kernel — 2026-09-18

Implementation `f77cd78`, base `1d637f3`.

The pinned PeleC boundary guide explicitly marks GC-NSCBC as not ported into
its C++ implementation. Its EB-ConvergingNozzle case does implement a
frozen-composition pressure outlet. This increment ports those equations,
generalized to all six Cartesian outward normals, into the public
`characteristic_pressure_ghost_3d` API in `simulation_config_reactive_3d_mod`.
It does not implement GC-NSCBC or connect a new solver boundary mode yet.

For ghost-distance ratio r (interior-to-ghost distance divided by
interior-to-boundary distance), delta=r*(target pressure-interior pressure).
Ghost density is rho+delta/c²; outward normal velocity is un-delta/(rho*c);
pressure is p+delta. Tangential velocities and mass fractions are preserved.
The first uniform-grid ghost center uses r=2; boundary-face state uses r=1.
NASA7 primitive conversion reconstructs temperature and total energy. Invalid
face/shape/finite/range inputs or nonphysical extrapolated states reject without
changing outputs. There is no automatic Mach classification or pressure clipping.

## Reference provenance

Fixed reference commit: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

- [Boundary guide](https://github.com/Pele-Suite/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Docs/sphinx/BoundaryConditions.rst), SHA256 `897b100e6ab931a850d0427dc324852ef3dfc98879772ba0c4c79351609c4f2d`.
- [Nozzle boundary implementation](https://github.com/Pele-Suite/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Exec/RegTests/EB-ConvergingNozzle/prob.H), outlet type 1, SHA256 `8e3d840a7857ec4d765019aad00369ea6fbe9aef5731fc3970c79c29e97c0bcb`.

Downloaded reference artifacts are in `build/nscbc-reference/`. No upstream
source was changed. Equations use interior EOS sound speed and frozen composition.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Test executed serially; selected
library compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_boundaries_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_boundaries_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, 0 failures/skips, 4.14 s total. Builds, provenance,
syntax, contract and diff checks passed. Six-face tests compare converted ghost
primitive values to the explicit formula at relative 3e-12 with 1e-30 floor,
including unchanged tangential velocities and species fractions. An extrapolated
negative pressure preserves ghost-state and temperature sentinels. Existing
PCM/PLM and three-Riemann-solver boundary tests remain.

Artifacts: `build/pressure-ghost-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
This is formula/kernel evidence, not a running C++ oracle, acoustic-reflection
qualification, boundary-flow validation or full characteristic inlet/outlet support.
