# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ syntax translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The current executable milestone provides:

- serial one-dimensional compressible Euler equations;
- constant-`gamma` ideal-gas EOS;
- conservative finite-volume updates;
- selectable Rusanov or qualified PeleC-style approximate Riemann fluxes;
- selectable `pcm`, componentwise `plm`, or time-traced `pelec_plm` reconstruction;
- minmod and monotonized-central (`mc`) limiting;
- outflow and periodic ghost-cell boundaries;
- SSPRK2 for `pcm`/`plm` and a single-stage time-centered Godunov update for `pelec_plm`;
- selectable Sod and Shu-Osher problems;
- smooth-convergence, exact-solution, deterministic-regression, conservation, and CI gates.

`pelec_plm` implements the one-dimensional, single-species, constant-`gamma` characteristic projection and wave tracing corresponding to the core structure of PeleC `Source/PLM.H`. It is intentionally qualified: fourth-order slopes, flattening, general-EOS terms, species tracing, embedded boundaries, and multidimensional transverse corrections are not implemented yet.

AMR, chemistry, diffusion, MPI, embedded boundaries, LES, and spray remain planned work.

## Build and test

Requirements: CMake 3.23 or newer, a Fortran 2018 compiler, and Python 3.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

With Ninja installed:

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

## Run representative cases

First-order Rusanov baseline:

```bash
./build/pelef cases/sod/sod.nml
python3 tools/compare_sod.py --input sod.csv
```

Componentwise PLM with the PeleC-style Riemann solver:

```bash
./build/pelef cases/sod/sod_pelec.nml
python3 tools/compare_sod.py --input sod_pelec.csv
```

PeleC-style characteristic tracing and Riemann solver:

```bash
./build/pelef cases/sod/sod_pelec_plm.nml
python3 tools/compare_sod.py \
  --input sod_pelec_plm.csv \
  --density-l1-max 1.6e-3 \
  --pressure-l1-max 1.0e-3
```

Shu-Osher with characteristic tracing:

```bash
./build/pelef cases/shu_osher/shu_osher_pelec_plm.nml
python3 tools/check_shu_osher.py \
  --input shu_osher_pelec_plm.csv \
  --density-squared-reference 112.78740201441508 \
  --density-moment-reference -7.071310814334509 \
  --signature-relative-tolerance 5e-6
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
