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

This older passive path intentionally retains the constant-`gamma` hydro baseline.

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

### `pelef_reactive_1d`: general-EOS reactive Euler solver

- conserved state `(rho, rho*u, rho*v, rho*w, rho*E, rho*Y_k)`;
- NASA7 composition-dependent pressure, temperature, heat capacities, ratio of specific heats, and frozen sound speed;
- safeguarded conserved-to-primitive recovery through `e(Y,T) -> T` inversion;
- PCM or frozen-composition characteristic PLM reconstruction;
- Rusanov flux with species-flux closure;
- periodic or outflow boundaries;
- cell-local adiabatic constant-volume chemistry;
- reaction-hydro-reaction Strang splitting;
- homogeneous-reactor reduction, smooth entropy-wave convergence, and nonuniform reactive-hotspot regressions.

The reactive path currently uses the verified seven-species, four-reaction elementary subset. The characteristic projection is a qualified frozen-composition ideal-gas-mixture approximation, not full PeleC/PelePhysics general-EOS characteristic parity.

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

Elementary H2/O2 constant-volume reactor:

```bash
./build/pelef0d_h2o2 cases/zero_d_h2o2/reactor.nml
python3 tools/check_zero_d_h2o2.py --input zero_d_h2o2.csv
```

Reactive one-dimensional hotspot:

```bash
./build/pelef_reactive_1d cases/reactive_hotspot/hotspot.nml
python3 tools/check_reactive_hotspot.py --input reactive_hotspot.csv
```

Smooth general-EOS entropy wave:

```bash
./build/pelef_reactive_1d \
  cases/reactive_entropy_wave/entropy_wave.nml
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
