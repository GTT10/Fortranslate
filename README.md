# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The `0.6.0` milestone contains two serial uniform-grid executables.

### `pelef`: one-dimensional Euler solver

- constant-`gamma` ideal-gas EOS;
- `pcm`, componentwise primitive `plm`, and time-traced characteristic `pelec_plm` reconstruction;
- order-2 or PeleC-style five-point order-4 limited slopes;
- optional pressure/velocity shock flattening;
- Rusanov or qualified single-species PeleC-style Riemann fluxes;
- outflow and periodic boundaries;
- Sod, Shu-Osher, and symmetric planar Sedov-type regressions.

### `pelef2d`: two-dimensional Euler scaffold

- uniform Cartesian mesh with periodic boundaries;
- shared five-component conserved state `(rho, rho*u, rho*v, rho*w, rho*E)`;
- x/y directional Riemann fluxes through an explicit momentum rotation;
- limited primitive slopes and normal characteristic tracing in both directions;
- conservative half-step transverse flux corrections;
- positivity-preserving scaling of transverse corrections when required;
- one unsplit conservative CTU-style update;
- periodic isentropic-vortex analytical and convergence regressions.

The two-dimensional implementation is a qualified regular-grid subset. It does not yet include PeleC multidimensional source terms, embedded boundaries, 3D double-transverse corrections, AMR, chemistry, transport, or species arrays.

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

One-dimensional characteristic PLM:

```bash
./build/pelef cases/sod/sod_pelec_plm.nml
python3 tools/compare_sod.py --input sod_pelec_plm.csv
```

Strong-blast regression with order-4 slopes and flattening:

```bash
./build/pelef cases/sedov/sedov.nml
python3 tools/check_sedov.py --input sedov.csv
```

Two-dimensional periodic isentropic vortex:

```bash
./build/pelef2d cases/isentropic_vortex/vortex.nml
python3 tools/check_isentropic_vortex.py \
  --input isentropic_vortex.csv
```

The 2D convergence test uses 24, 48, and 96 cells per direction. The observed density orders are approximately `2.278` and `2.276`. At 48 cells, enabling the transverse correction reduces density L1 error from approximately `1.149e-3` to `5.333e-4`.

## Project records

- [Porting plan](docs/porting_plan.md)
- [Architecture](docs/architecture.md)
- [Numerical methods](docs/numerical_methods.md)
- [PeleC responsibility mapping](docs/pelec_mapping.md)
- [State variables](docs/state_variables.md)
- [Parity strategy](docs/parity_strategy.md)
- [Implementation status](docs/implementation_status.md)
- [Design decisions](docs/design_decisions/)
