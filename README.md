# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The `0.52.0` milestone contains ten serial verification executables, six
optional MPI verification executables, and a runnable one-dimensional reactive
AMR application with solution-driven dynamic regridding and molecular
transport.

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

The reactive applications can select either the verified seven-species,
four-reaction elementary subset or the full ten-species, 29-reaction H2/O2
mechanism with third-body, falloff, Troe, and adaptive implicit chemistry. The
characteristic projection is a qualified frozen-composition ideal-gas-mixture
approximation, not full PeleC/PelePhysics general-EOS characteristic parity.

### MPI one-dimensional verification

With `PELEF_ENABLE_MPI=ON`, five existing executables verify:

- uneven non-replicated block decomposition and periodic halo exchange;
- conservative multispecies Euler transport;
- distributed adaptive implicit full-H2/O2 chemistry scheduling;
- general-EOS molecular transport with viscosity, conduction, barodiffusion,
  correction velocity, and species-enthalpy transport;
- transactional reaction--transport--hydro--transport--reaction splitting;
- ordered gather output, global timestep/conservation reductions, and
  complete-field parity for 1, 2, 4, and 8 ranks.

A sixth MPI executable exercises the AMR distribution bridge. It validates a
replicated patch-tree description collectively, assigns every root/fine patch
to one deterministic cell-weighted owner, broadcasts owner-authoritative
patch fields, and exchanges up to four adjacent-sibling halo layers across
rank boundaries. It now also advances chemistry on owners alone, reaches
communicator-wide success/failure consensus after every patch, broadcasts the
accepted reactive state, averages down deepest-to-root, rebuilds ghosts, and
rolls every rank back exactly after a rejected owner update. It also drives
recursive hyperbolic advancement on patch owners alone, broadcasts each
owner's start state and complete face-flux field, then reuses the serial child
subcycling, fine/fine flux reconciliation, reflux, average-down, and ghost
rules on synchronized replicas. The current bridge retains replicated field
storage. Molecular transport now follows the same owner-only recursion with
parabolic `r^2` subcycling and shared diffusive fluxes. A single distributed
full-physics transaction and sparse rank-local patch storage remain the next
integration steps.

### One-dimensional AMR

The AMR layer provides:

- a static two-level hierarchy with an interior refined patch;
- an arbitrary-depth hierarchy foundation built from nested level interfaces;
- arbitrary component counts and integer refinement ratios;
- MC-limited, conservative piecewise-linear prolongation;
- volume-average restriction and covered-cell synchronization;
- refinement-ratio level subcycling;
- time-integrated coarse/fine flux registers and reflux;
- a composite-integral gate proving conservation across both patch interfaces;
- component-selectable normalized-gradient tagging with absolute noise floors;
- buffered, minimum-width single-patch planning with explicit boundary checks;
- conservative patch creation, movement, resizing, and removal;
- exact retention of same-resolution fine data where old and new patches overlap;
- reactive general-EOS states advanced on both coarse and fine levels;
- refinement-ratio fine-level subcycling with time-interpolated coarse ghosts;
- coarse/fine flux accumulation, reflux, and average-down every coarse step;
- symmetric chemistry--hydro--chemistry splitting across the AMR hierarchy;
- transactional hierarchy rollback and periodic solution-driven regridding;
- ordered composite CSV output with exact domain-coverage checks;
- optional limited primitive-variable PLM with SSPRK2 time advancement;
- time-averaged SSPRK2 interface fluxes used consistently for reflux;
- four-layer physical and coarse/fine ghost storage for PPM stencils;
- primitive or time-traced characteristic PPM with SSPRK3 level advancement;
- optional PeleC-style hybrid WENO5-JS, WENO5-Z, WENO7-Z, or WENO3-Z edge
  reconstruction inside the characteristic-PPM path;
- SSPRK3 effective interface fluxes used consistently for reflux;
- conservative, MC-limited parent-to-fine ghost interpolation at subcycle
  midpoint times;
- AMR viscosity, Fourier conduction, and mixture-averaged species diffusion;
- parabolic fine subcycling with time-interpolated coarse transport ghosts;
- diffusive flux-register reflux and covered-cell average-down;
- symmetric reaction--transport--hydro--transport--reaction composition;
- cumulative multilevel subcycle schedules and deepest-to-root synchronization;
- a four-level, mixed-ratio conservation gate across every interface;
- arbitrary-depth reactive state and temperature ownership;
- recursive hydro subcycling and diffusive `r^2` subcycling at every depth;
- recursive chemistry and reaction--transport--hydro splitting;
- a three-level reactive conservation, positivity, closure, and synchronization
  gate;
- a runtime `amr_max_levels` limit with tag-driven nested hierarchy creation;
- periodic conservative multilevel hierarchy rebuilds during simulation;
- recursively ordered, exact-coverage multilevel composite CSV output;
- a runnable three-level hotspot case and structural output gate;
- physical-coordinate overlap transfer across changed multilevel hierarchies;
- exact retention of aligned old fine state and temperature data;
- nested fine patches touching an outflow physical boundary;
- physical-side PPM/WENO ghost fill with reflux restricted to the remaining
  coarse/fine interface;
- ordered sets of disjoint fine patches over one parent level;
- disconnected-tag clustering with deterministic buffer/minimum-width
  expansion and automatic coalescing of adjacent candidates;
- set-wide conservative prolongation, average-down, reflux, and composite
  integration without double counting covered parent cells;
- conservative patch-set movement, repartition, removal, and exact retention
  of same-resolution fine overlap across old/new patch intersections;
- fixed two-level reactive WENO7-Z hydro subcycling on two separated fine
  patches with per-patch flux registers;
- fixed two-level multipatch reaction--transport--hydro splitting with `r^2`
  diffusive subcycling, per-patch reflux, and transactional rollback;
- tag-driven two-level multipatch creation, movement, repartition, removal,
  overlap retention, runtime statistics, and ordered composite CSV output;
- arbitrary-depth separated patch trees with explicit parent ownership,
  mixed per-level refinement ratios, recursive conservative prolongation,
  deepest-to-root average-down, and exact composite integration;
- parent-owned patch-tree flux registers with transactional deepest-to-root
  reflux and covered-cell synchronization across every branch;
- static arbitrary-depth reactive patch trees with recursive PCM hydro
  subcycling, time-interpolated parent ghosts, per-child flux accumulation,
  reflux, and average-down at every branch;
- a four-level `1/2/3/2`-patch reactive conservation, synchronization,
  positivity, closure, and exact subcycle-count gate;
- symmetric chemistry--hydro--chemistry splitting over every patch-tree node,
  with deepest-to-root post-reaction synchronization and whole-tree rollback;
- a chemistry-on versus hydro-only branch comparison proving reactive species
  evolution while retaining composite mass, momentum, and total energy;
- recursive molecular transport with cumulative `r^2` child subcycling,
  time-interpolated parent ghosts, per-child diffusive flux registers, reflux,
  and average-down at every branch;
- full reaction--transport--hydro--transport--reaction composition with a
  four-level transport-on versus transport-off conservation gate;
- transactional runtime patch-tree rebuilds from supplied branching plans,
  with physical-coordinate overlap transfer across changed parent ownership;
- no-op, moved-tree conservation, deepest exact state/temperature retention,
  counter preservation, and invalid-plan rollback gates;
- per-parent normalized-gradient tagging and deterministic disconnected-tag
  clustering for automatic arbitrary-depth branching plans;
- tag-driven transactional tree rebuilds, including root-only creation,
  maximum-depth branching, unchanged-plan no-op, and invalid-request gates;
- independently owned adjacent children with parent-local same-level exchange
  for face and four-layer PPM/WENO ghost data;
- single owned time-integrated hydro and diffusive fluxes at fine/fine faces,
  with those internal faces excluded from coarse/fine reflux;
- adjacent PPM hydro and molecular-transport conservation, synchronization,
  exact exchange, and subcycle-accounting gates;
- deterministic cell-weighted MPI ownership for every tree patch, collective
  hierarchy-consensus rejection, owner-authoritative patch synchronization,
  and four-layer cross-rank adjacent-sibling halo gates;
- owner-only patch-tree chemistry with one global advance per patch, serial
  field parity, deepest-to-root synchronization, conservation, and global
  transactional rollback gates;
- owner-only recursive MPI patch-tree hydro with exact per-owner subcycle
  accounting, cross-owner adjacent PPM flux reconciliation, serial field
  parity, conservation, and global transactional rollback gates;
- owner-only recursive MPI patch-tree molecular transport with exact
  parabolic subcycle accounting, cross-owner shared diffusive fluxes, serial
  parity, conservation, and global transactional rollback gates;
- a moving-contact gate demonstrating lower AMR error than PCM.

For PCM/PLM, the reactive AMR application retains its overlap-preserving
two-level path when `amr_max_levels = 2`. PPM selects the multilevel engine at
any configured depth so every level owns the wider stencil state. The same
engine provides tag-driven arbitrary-depth state ownership, recursive
advancement, periodic hierarchy rebuilds, and composite output. A changed
multilevel hierarchy is conservatively averaged to the root before nested
patches are rebuilt, then old fine state and temperature data are copied
wherever old/new level spacing and physical cells align. Changed refinement
ratios fall back to conservative prolongation. With
`amr_multipatch_enabled = .true.`, the public application uses a tag-driven
two-level patch set, clusters disconnected tags, periodically rebuilds the
set, retains aligned fine overlap, and writes every uncovered parent or fine
cell exactly once. A separate static patch-tree engine permits each parent
patch to own zero or more separated children at arbitrary depth and advances
those patches recursively with level-ratio hydro subcycling. Every local child
accumulates its own coarse/fine boundary fluxes; reflux and average-down occur
after its fine subcycles, and failures restore the complete tree. The
patch-tree engine can synchronize to the root, tag and cluster every
prospective parent independently, and rebuild the resulting arbitrary-depth
branching plan transactionally. Independently owned adjacent siblings exchange
same-level ghosts, reconcile each shared time-integrated interface flux, and
exclude that internal side from reflux. The first MPI distribution bridge
assigns a unique owner to every patch and communicates authoritative fields
and adjacent halos while retaining replicas on all ranks. Chemistry and the
recursive hydro patch kernel execute only on those owners with collective
acceptance and rollback. Owner-authoritative start-state and face-flux
broadcasts let every replica apply the existing subcycling, shared-flux,
reflux, and average-down rules deterministically. Owner-only molecular
transport uses the same structure with `r^2` subcycling. A single distributed
full-physics transaction and sparse patch storage remain later slices. A
periodic child may
touch a physical boundary only when it covers the full parent domain;
one-sided periodic refinement remains excluded because it crosses the periodic
seam.

## Build and test

Requirements: CMake 3.23 or newer, a Fortran 2018 compiler, and Python 3.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

MPI verification additionally requires an MPI implementation with the Fortran
2018 `mpi_f08` module:

```bash
cmake -S . -B build-mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_MPI=ON
cmake --build build-mpi --parallel
mpiexec -n 4 ./build-mpi/pelef_mpi_reactive_1d mpi_reactive_np4.csv
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


### PeleF 0.24.0 MPI 1D verification

Configure with `-DPELEF_ENABLE_MPI=ON`, then run the MPI verification drivers
with 1, 2, 4, or 8 ranks. The foundation and multispecies drivers use 257 cells
so the block decomposition is intentionally uneven; the chemistry, transport,
and coupled-reactive drivers use smaller non-divisible workloads to exercise the
same decomposition and ordered-gather logic.
