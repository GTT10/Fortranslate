# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ syntax translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The current executable milestone provides:

- serial one-dimensional compressible Euler equations;
- constant-gamma ideal-gas EOS;
- conservative finite-volume updates;
- selectable Rusanov or single-species PeleC-style approximate Riemann fluxes;
- selectable piecewise-constant (`pcm`) or piecewise-linear (`plm`) reconstruction;
- selectable minmod or monotonized-central (`mc`) slope limiting;
- outflow and periodic ghost-cell boundaries;
- SSPRK2 time integration;
- selectable Sod and Shu-Osher problems;
- smooth periodic entropy-wave convergence tests;
- unit, exact-solution, deterministic-regression, conservation, and CI gates.

The PLM implementation is componentwise in primitive variables. The PeleC-style Riemann solver is a constant-`gamma`, single-species reduction of the acoustic star-state and wave-interpolation logic in `Source/Riemann.H`. These are verified intermediate components, not claims of complete PeleC Godunov parity.

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

## Run the hydro cases

First-order Rusanov baseline:

```bash
./build/pelef cases/sod/sod.nml
python3 tools/compare_sod.py --input sod.csv
```

PLM/MC with Rusanov:

```bash
./build/pelef cases/sod/sod_plm.nml
python3 tools/compare_sod.py \
  --input sod_plm.csv \
  --density-l1-max 4e-3 \
  --pressure-l1-max 3e-3
```

PLM/MC with the PeleC-style approximate Riemann solver:

```bash
./build/pelef cases/sod/sod_pelec.nml
python3 tools/compare_sod.py \
  --input sod_pelec.csv \
  --density-l1-max 2e-3 \
  --pressure-l1-max 1.5e-3
```

Shu-Osher shock-density-wave interaction:

```bash
./build/pelef cases/shu_osher/shu_osher.nml
python3 tools/check_shu_osher.py --input shu_osher.csv
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
