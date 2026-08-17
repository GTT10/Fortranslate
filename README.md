# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The current `0.5.0` milestone provides a serial one-dimensional solver for the compressible Euler equations with a constant-`gamma` ideal-gas EOS.

Available numerical components:

- reconstruction: `pcm`, componentwise primitive `plm`, or time-traced characteristic `pelec_plm`;
- PLM slope order for `pelec_plm`: `plm_order = 2` or `4`;
- optional PeleC pressure/velocity shock flattening with `use_flattening = .true.`;
- slope limiter: `minmod` or monotonized-central (`mc`);
- Riemann solver: `rusanov` or the documented single-species PeleC-style subset;
- boundaries: `outflow` or `periodic`;
- problems: Sod, Shu-Osher, and a symmetric planar Sedov-type blast.

The fourth-order slope and flattening formulas follow the corresponding one-dimensional regular-cell logic in PeleC `Source/PLM.H` and `Source/Godunov.H`. General-EOS terms, species tracing, embedded boundaries, and multidimensional transverse corrections are not yet implemented.

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

## Representative runs

Characteristic PLM with the qualified PeleC-style Riemann solver:

```bash
./build/pelef cases/sod/sod_pelec_plm.nml
python3 tools/compare_sod.py --input sod_pelec_plm.csv
```

Fourth-order characteristic PLM with flattening on the strong-blast case:

```bash
./build/pelef cases/sedov/sedov.nml
python3 tools/check_sedov.py --input sedov.csv
```

The Sedov regression checks positivity, reflection symmetry, integral conservation, shock location, and pinned field signatures.

## Project records

- [Porting plan](docs/porting_plan.md)
- [Architecture](docs/architecture.md)
- [Numerical methods](docs/numerical_methods.md)
- [PeleC responsibility mapping](docs/pelec_mapping.md)
- [State variables](docs/state_variables.md)
- [Parity strategy](docs/parity_strategy.md)
- [Implementation status](docs/implementation_status.md)
- [Design decisions](docs/design_decisions/)
