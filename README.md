# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ syntax translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The current executable milestone provides:

- serial one-dimensional compressible Euler equations;
- constant-gamma ideal-gas EOS;
- conservative finite-volume updates;
- Rusanov/local Lax–Friedrichs interface fluxes;
- selectable piecewise-constant (`pcm`) or piecewise-linear (`plm`) reconstruction;
- selectable minmod or monotonized-central (`mc`) slope limiting;
- outflow and periodic ghost-cell boundaries;
- SSPRK2 time integration;
- canonical Sod shock-tube cases;
- a smooth periodic entropy-wave convergence test;
- unit tests, exact-solution comparisons, conservation checks, and CI.

The PLM implementation is a verified componentwise primitive-variable scheme. It is an intermediate step toward PeleC-style characteristic reconstruction and tracing; it is not yet a claim of full `Source/PLM.H` parity.

AMR, chemistry, diffusion, MPI, embedded boundaries, LES, and spray are planned but are not implemented yet.

## Build and test

Requirements: CMake 3.23 or newer, a Fortran 2018 compiler, and Python 3 for regression comparisons.

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

## Run the Sod cases

First-order baseline:

```bash
./build/pelef cases/sod/sod.nml
python3 tools/compare_sod.py --input sod.csv
```

PLM with the MC limiter:

```bash
./build/pelef cases/sod/sod_plm.nml
python3 tools/compare_sod.py \
  --input sod_plm.csv \
  --density-l1-max 4e-3 \
  --pressure-l1-max 3e-3
```

When using the debug preset, the executable path is `build/debug/pelef`.

## Project records

- [Porting plan](docs/porting_plan.md)
- [Architecture](docs/architecture.md)
- [Numerical methods](docs/numerical_methods.md)
- [PeleC responsibility mapping](docs/pelec_mapping.md)
- [State variables](docs/state_variables.md)
- [Parity strategy](docs/parity_strategy.md)
- [Implementation status](docs/implementation_status.md)
- [Design decisions](docs/design_decisions/)
