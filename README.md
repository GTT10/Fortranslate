# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The `0.19.0` milestone contains eight serial verification executables.

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


### `pelef_transport_probe`: dilute-gas transport verification

- Lennard-Jones transport records for H2, H, O, O2, OH, H2O, and N2;
- Chapman--Enskog pure viscosities and binary diffusion coefficients;
- Wilke mixture viscosity;
- modified-Eucken pure conductivity and Mathur mixture conductivity;
- mixture-averaged species diffusion coefficients;
- live qualification against Cantera 3.2 at four temperature, pressure, and
  composition states.

This is a deliberately qualified ideal-gas subset. It does not yet reproduce
PelePhysics polynomial transport fits, polar corrections, Soret diffusion,
multicomponent diffusion, or bulk viscosity.

### `pelef_reactive_1d`: general-EOS reactive Euler solver

- conserved state `(rho, rho*u, rho*v, rho*w, rho*E, rho*Y_k)`;
- NASA7 composition-dependent pressure, temperature, heat capacities, ratio of specific heats, and frozen sound speed;
- safeguarded conserved-to-primitive recovery through `e(Y,T) -> T` inversion;
- PCM, frozen-composition characteristic PLM, monotone primitive PPM, or
  time-traced frozen-composition characteristic PPM reconstruction;
- selectable Rusanov or general-EOS HLLC flux with species-flux closure;
- SSPRK3 time integration for the semidiscrete primitive-PPM path;
- PeleC-style parabolic profile integration over the `u-c`, `u`, and `u+c`
  waves for characteristic PPM;
- optional PeleC one-dimensional shock flattening and bounded
  Colella--Woodward contact steepening on characteristic PPM;
- periodic or outflow boundaries;
- cell-local adiabatic constant-volume chemistry;
- optional shear viscosity, Fourier heat conduction, mixture-averaged species
  diffusion, barodiffusion, correction velocity, and species-enthalpy flux;
- explicit SSPRK2 transport with a parabolic timestep gate;
- symmetric reaction--transport--hydro--transport--reaction splitting;
- homogeneous-reactor reduction, smooth density/composition-wave convergence,
  discontinuous material-contact resolution, primitive/characteristic PPM
  convergence, periodic strong-shock flattening, and nonuniform
  reactive-hotspot regressions.


### `pelef_reactive_2d`: general-EOS reactive CTU solver

- uniform periodic Cartesian mesh with the same conserved reactive state as the 1D path;
- composition-dependent NASA7 pressure, temperature, heat capacities, and frozen sound speed;
- directional general-EOS Rusanov or HLLC fluxes through explicit momentum rotation;
- PCM, frozen-composition characteristic PLM, or time-traced
  frozen-composition characteristic PPM in both coordinate directions;
- optional bounded contact steepening and PeleC-style shock flattening on the
  characteristic-PPM normal predictor;
- provisional face fluxes, conservative CTU transverse half-step corrections, and EOS-based positivity scaling;
- species, momentum, and total energy corrected together so `sum(rho*Y_k)=rho` remains coupled to the hydro update;
- cell-local chemistry with symmetric reaction--transport--hydro--transport--reaction splitting;
- optional x/y Newtonian viscosity, Fourier conduction, mixture-averaged species diffusion, barodiffusion, correction velocity, and species-enthalpy flux;
- exact diagonal density/composition-wave convergence, x/y one-dimensional
  reduction, material-contact sharpening, oblique strong-shock flattening,
  periodic vortex, and reacting-hotspot regressions.

The reactive 2D path supports matched periodic pairs, slip/no-slip walls,
adiabatic/isothermal wall temperatures, fixed-state inflow, and zero-gradient
outflow. Periodic cases retain the qualified characteristic-PPM plus CTU path;
physical faces use boundary-aware ghost reconstruction, an exact impermeable
inviscid wall flux, and boundary-aware molecular transport.

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

Reactive one-dimensional hotspot with characteristic PLM:

```bash
./build/pelef_reactive_1d cases/reactive_hotspot/hotspot.nml
python3 tools/check_reactive_hotspot.py --input reactive_hotspot.csv
```

The same case with monotone PPM and HLLC:

```bash
./build/pelef_reactive_1d cases/reactive_hotspot/hotspot_ppm.nml
python3 tools/check_reactive_hotspot.py --input reactive_hotspot_ppm.csv
```

The time-traced characteristic PPM path with the optional contact and shock
detectors enabled:

```bash
./build/pelef_reactive_1d \
  cases/reactive_hotspot/hotspot_characteristic_ppm.nml
python3 tools/check_reactive_hotspot.py \
  --input reactive_hotspot_characteristic_ppm.csv
```

Smooth general-EOS entropy wave:

```bash
./build/pelef_reactive_1d \
  cases/reactive_entropy_wave/entropy_wave.nml
```

General-EOS H2/N2 composition wave with HLLC:

```bash
./build/pelef_reactive_1d \
  cases/reactive_composition_wave/composition_wave.nml
python3 tools/check_reactive_composition_wave.py \
  --input reactive_composition_wave.csv
```


Dilute-gas transport coefficient probe:

```bash
./build/pelef_transport_probe transport_probe.csv
python3 tools/compare_transport_cantera.py --input transport_probe.csv
```

Periodic one-dimensional molecular-transport pulse:

```bash
./build/pelef_reactive_1d \
  cases/reactive_transport_1d/transport_pulse.nml
python3 tools/check_reactive_transport_1d.py \
  --input reactive_transport_pulse.csv --nx 96
```

Reactive two-dimensional hotspot with CTU and HLLC:

```bash
./build/pelef_reactive_2d cases/reactive_hotspot_2d/hotspot.nml
python3 tools/check_reactive_hotspot_2d.py \
  --input reactive_hotspot_2d.csv --nx 24 --ny 24
```

The same two-dimensional hotspot with characteristic PPM and the optional
contact/shock controls enabled:

```bash
./build/pelef_reactive_2d \
  cases/reactive_hotspot_2d/hotspot_characteristic_ppm.nml
python3 tools/check_reactive_hotspot_2d.py \
  --input reactive_hotspot_characteristic_ppm_2d.csv --nx 24 --ny 24
```

Oblique exact entropy-wave transport through the same 2D path:

```bash
./build/pelef_reactive_2d \
  cases/reactive_diagonal_wave_2d/diagonal_wave.nml
```

Oblique constant-pressure H2/N2 composition transport through the
characteristic-PPM/CTU path:

```bash
./build/pelef_reactive_2d \
  cases/reactive_diagonal_wave_2d/diagonal_composition_ppm.nml
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


### PeleF 0.17.0 two-dimensional transport example

```bash
./build/pelef_reactive_2d cases/reactive_transport_2d/transport_hotspot.nml
python3 tools/check_reactive_hotspot_2d.py \
  --input reactive_transport_hotspot_2d.csv --nx 20 --ny 20
```


### PeleF 0.18.0 physical-boundary examples

```bash
./build/pelef_reactive_2d cases/reactive_boundaries_2d/couette.nml
./build/pelef_reactive_2d cases/reactive_boundaries_2d/thermal_channel.nml
./build/pelef_reactive_2d cases/reactive_boundaries_2d/inflow_outflow.nml
```

Solid walls are species-impermeable. Slip walls remove tangential viscous
stress; no-slip walls reflect velocity about a prescribed wall velocity.


### Full pressure-dependent H2/O2 chemistry

Set `chemistry_model = "full_h2o2"` to use the 10-species, 29-reaction third-body/Troe mechanism with the implicit cell reactor.


### PeleF 0.20.0 MPI 1D verification

Configure with `-DPELEF_ENABLE_MPI=ON`, then run `pelef_mpi_1d` with 1, 2, or 4 ranks. The verification driver uses 257 cells so the block decomposition is intentionally uneven.
