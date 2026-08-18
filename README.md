# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The `0.8.0` milestone contains four serial executables.

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
- x/y directional Riemann fluxes through explicit momentum rotation;
- limited primitive slopes and normal characteristic tracing;
- CTU-style transverse half-step corrections with positivity scaling;
- one unsplit conservative update;
- periodic isentropic-vortex analytical and convergence regressions.

### `pelef_ms`: passive multispecies Euler solver

- runtime species count with conserved `rho*Y_k` components;
- checked mass-fraction conversion, positivity, and species closure;
- species fluxes tied exactly to the shared mass flux;
- one-dimensional characteristic tracing and two-dimensional CTU transport;
- MultiSpecSod and periodic species-wave regressions.

This hydro path remains passive and constant-`gamma`: composition does not yet change the flux-level pressure, heat capacity, sound speed, or temperature.

### `pelef0d`: thermodynamics and reactor verification

- species molecular weights and NASA7 thermodynamic polynomials;
- mass-based mixture molecular weight, gas constant, `cp`, `cv`, `gamma`, enthalpy, internal energy, and frozen sound speed;
- ideal-gas pressure/density conversion;
- bracketed Newton/bisection inversion from specific internal energy to temperature;
- a constant-volume two-species isomerization reactor;
- isothermal analytical and adiabatic energy-conservation gates.

The built-in H2, O2, H2O, and N2 coefficients are a small verified subset of the Cantera GRI-Mech/air data. The reactor mechanism is deliberately synthetic and is not a detailed combustion model.

AMR, embedded boundaries, diffusion, detailed chemistry, mechanism parsing, MPI, accelerator execution, LES, and spray remain future work.

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

Two-dimensional periodic isentropic vortex:

```bash
./build/pelef2d cases/isentropic_vortex/vortex.nml
python3 tools/check_isentropic_vortex.py --input isentropic_vortex.csv
```

Passive two-species Sod problem:

```bash
./build/pelef_ms cases/multispec_sod/multispec_sod.nml
python3 tools/check_multispec_sod.py --input multispec_sod.csv
```

Adiabatic zero-dimensional isomerization reactor:

```bash
./build/pelef0d cases/zero_d_isomerization/reactor.nml
python3 tools/check_zero_d_isomerization.py \
  --input zero_d_isomerization.csv
```

## Project records

- [Porting plan](docs/porting_plan.md)
- [Architecture](docs/architecture.md)
- [Numerical methods](docs/numerical_methods.md)
- [PeleC responsibility mapping](docs/pelec_mapping.md)
- [State variables](docs/state_variables.md)
- [Parity strategy](docs/parity_strategy.md)
- [Implementation status](docs/implementation_status.md)
- [Design decisions](docs/design_decisions/)
