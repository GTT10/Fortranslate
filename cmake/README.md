# CMake layout

The root `CMakeLists.txt` owns the project version, options, dependency checks,
and compiler standard. Its remaining includes intentionally preserve the
original registration order and directory scope:

1. `targets/Core.cmake`: fixed-mechanism core and compiler options.
2. `targets/SelectedRuntime.cmake`: isolated selected-mechanism runtimes and
   the optional official-Fortran CVODE backend.
3. `targets/Applications.cmake`: selected/fixed serial programs and installation.
4. `add_subdirectory(tests)`: test registration, including selected fixtures.
5. `targets/MPI.cmake`: fixed MPI support/programs and their tests/installation.
6. `targets/InstalledTests.cmake`: private-prefix installed restart tests.

`tests/CMakeLists.txt` includes 19 topic-specific files under `tests/cmake/`.
Their order is part of the dependency contract. Some later registrations use
fixtures, directories, or test dependencies introduced by earlier files.
Do not sort the includes alphabetically or replace `include()` with
`add_subdirectory()` without checking the changed directory/variable scope.

Keep new registrations in the corresponding module. Put reusable CMake
functions beside `MechanismBundle.cmake`; do not restore a monolithic root file.
The split introduces no new numerical functionality, compiler flags, runtime
options, test tolerances, expected hashes, or installation destinations.

The migration evidence is in
[the layout validation record](../docs/validation/build-layout-20260909.md).
It compares generated target graphs, CTest commands/properties, and install
rules across seven configurations. Numerical CI remains a separate gate.
