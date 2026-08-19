# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The `0.10.0` milestone contains six serial verification executables.

### `pelef`: one-dimensional Euler solver

- constant-`gamma` ideal-gas EOS;
- PCM, componentwise primitive PLM, and time-traced characteristic PLM;
- order-2 or PeleC-style five-point order-4 limited slopes;
- optional pressure/velocity shock flattening;
- Rusanov or qualified single-species PeleC-style Riemann fluxes;
- outflow and periodic boundaries;
- Sod, Shu-Osher, and symmetric planar Sedov-type regressions.

### `pelef2d`: two-dimensional Euler scaffold

- uniform periodic Cartesian mesh;
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

This hydro path remains passive and constant-`gamma`: composition does not yet alter flux-level pressure, heat capacity, sound speed, or temperature.

### `pelef0d`: thermodynamics and toy-reactor verification

- species molecular weights and NASA7 thermodynamic polynomials;
- mass-based mixture molecular weight, gas constant, `cp`, `cv`, `gamma`, enthalpy, internal energy, and frozen sound speed;
- ideal-gas pressure/density conversion;
- bracketed Newton/bisection inversion from specific internal energy to temperature;
- a synthetic constant-volume two-species isomerization reactor;
- isothermal analytical and adiabatic energy-conservation gates.

### `pelef0d_h2o2`: generated elementary H2/O2 kinetics

- runtime elementary-reaction records with arbitrary reactant/product stoichiometry;
- reversible Arrhenius rates and NASA7 equilibrium constants;
- molar concentrations, progress rates, production rates, and mass-fraction source terms;
- JSON-to-Fortran mechanism generation with a committed-source cleanliness gate;
- an adaptive explicit RK4 constant-volume, adiabatic reactor;
- a seven-species, four-reaction H2/O2/N2 subset selected from Cantera `h2o2.yaml`;
- live trajectory and exact-state production-rate comparison against Cantera 3.2.

The seven-species case remains a deliberately limited elementary subset and is retained as a fast, high-precision kinetics regression.

### `pelef0d_h2o2_full`: pressure-dependent full H2/O2 verification

- ten species and all 29 reactions from Cantera `h2o2.yaml`;
- elementary, third-body, duplicate, and Troe falloff reaction forms;
- species-specific third-body efficiencies;
- analytic concentration and mass-fraction Jacobian assembly;
- reduced constant-energy reactor Jacobian including the temperature-composition coupling;
- safeguarded backward Euler/Newton solves with adaptive step doubling;
- Richardson extrapolation of accepted steps;
- structural conservation checks and live Cantera 3.2 trajectory/rate parity.

The implicit solver is a dense in-tree verification integrator, not CVODE. Direct Cantera/CHEMKIN parsing, sparse linear algebra, hydrocarbon mechanisms, and chemistry-flow coupling remain future work.

## Build and test

Requirements: CMake 3.23 or newer, a Fortran 2018 compiler, and Python 3.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

To enable the live Cantera reference gate:

```bash
python3 -m pip install cantera==3.2.0
cmake -S . -B build-cantera \
  -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_CANTERA_REFERENCE=ON
cmake --build build-cantera --parallel
ctest --test-dir build-cantera --output-on-failure
```

With Ninja installed, the provided presets can be used for the ordinary suite:

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

Synthetic adiabatic isomerization reactor:

```bash
./build/pelef0d cases/zero_d_isomerization/reactor.nml
python3 tools/check_zero_d_isomerization.py \
  --input zero_d_isomerization.csv
```

Elementary H2/O2 constant-volume reactor:

```bash
./build/pelef0d_h2o2 cases/zero_d_h2o2/reactor.nml
python3 tools/check_zero_d_h2o2.py --input zero_d_h2o2.csv
```

With Cantera installed, the same CSV can be compared directly:

```bash
python3 tools/compare_h2o2_cantera.py \
  --input zero_d_h2o2.csv \
  --mechanism mechanisms/h2o2_elementary_cantera.yaml
```

Full pressure-dependent H2/O2 reactor:

```bash
./build/pelef0d_h2o2_full cases/zero_d_h2o2_full/reactor.nml
python3 tools/check_zero_d_h2o2_full.py \
  --input zero_d_h2o2_full.csv
```

The live full-mechanism Cantera comparison is enabled through `PELEF_ENABLE_CANTERA_REFERENCE`.

## Project records

- [Porting plan](docs/porting_plan.md)
- [Architecture](docs/architecture.md)
- [Numerical methods](docs/numerical_methods.md)
- [PeleC responsibility mapping](docs/pelec_mapping.md)
- [State variables](docs/state_variables.md)
- [Parity strategy](docs/parity_strategy.md)
- [Implementation status](docs/implementation_status.md)
- [Design decisions](docs/design_decisions/)
