# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ syntax translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The repository now contains the first executable milestone:

- serial one-dimensional compressible Euler solver;
- constant-gamma ideal-gas EOS;
- conservative finite-volume update;
- Rusanov/local Lax–Friedrichs interface flux;
- SSPRK2 time integration;
- canonical Sod shock-tube case;
- unit tests, exact-solution comparison, conservation checks, and CI.

AMR, chemistry, diffusion, MPI, embedded boundaries, LES, and spray are planned but are not implemented yet.

## Build and test

Requirements: CMake 3.23 or newer, a Fortran 2018 compiler, and Python 3 for the regression comparison.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

With Ninja installed, the provided presets can be used:

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

## Run the Sod case

```bash
./build/pelef cases/sod/sod.nml
python3 tools/compare_sod.py --input sod.csv
```

When using the debug preset, the executable path is `build/debug/pelef`.

## Project records

- [Porting plan](docs/porting_plan.md)
- [Architecture](docs/architecture.md)
- [PeleC responsibility mapping](docs/pelec_mapping.md)
- [State variables](docs/state_variables.md)
- [Parity strategy](docs/parity_strategy.md)
- [Implementation status](docs/implementation_status.md)

Active development branch: `agent/pelec-fortran-port`.
