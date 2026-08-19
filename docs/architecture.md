# PeleF architecture

## Executable split

PeleF exposes six serial verification drivers over shared numerical and physical-property modules.

```text
pelef
  └─ one-dimensional constant-gamma PCM / PLM / PeleC-style paths

pelef2d
  └─ two-dimensional periodic constant-gamma CTU-style Euler path

pelef_ms
  └─ passive runtime-multispecies transport over the constant-gamma core

pelef0d
  └─ NASA7 mixture thermodynamics and a synthetic isomerization reactor

pelef0d_h2o2
  └─ generated reversible elementary kinetics and an H2/O2 reactor

pelef_reactive_1d
  └─ NASA7 general-EOS reactive Euler with characteristic PLM and Strang splitting
```

The established constant-`gamma` solvers remain intact as regression baselines. Composition-dependent flow is introduced through a separate driver and modules so a general-EOS change cannot silently alter earlier results.

## Shared thermodynamics and chemistry

```text
thermo_database_mod
  └─ named NASA7 species records and molecular weights
        ↓
nasa7_thermo_mod
  └─ species cp, cv, h, u, and standard-state entropy
        ↓
mixture_thermo_mod
  ├─ mixture molecular weight and gas constant
  ├─ mixture cp, cv, gamma, h, and u
  ├─ ideal-gas pressure, density, and frozen sound speed
  └─ safeguarded e(Y,T) -> T inversion

mechanisms/h2o2_elementary.json
        ↓
tools/generate_elementary_mechanism.py
        ↓
src/generated/h2o2_elementary_mechanism_mod.F90
        ↓
elementary_kinetics_mod
        ↓
constant_volume_reactor_mod
```

The normalized JSON file is the authoring format for the current generated reaction subset. CI regenerates the committed Fortran module and requires byte-for-byte equality.

## Reactive one-dimensional path

The reactive state is stored as

```text
state(variable, 0:nx+1)
```

with the conserved layout

```text
rho, rho*u, rho*v, rho*w, rho*E, rho*Y_1 ... rho*Y_N.
```

Temperature is a synchronized derived field, not an independently evolved conserved variable.

```text
state + temperature guess
        ↓
reactive_conserved_to_primitive
  ├─ recover Y from rho*Y
  ├─ remove kinetic energy from rho*E
  ├─ solve u(Y,T) = e_target
  ├─ evaluate p(Y,rho,T)
  └─ evaluate frozen sound speed
        ↓
reconstruction / Riemann / CFL
```

Hydrodynamic responsibilities are separated as follows:

```text
reactive_primitive_to_conserved
reactive_conserved_to_primitive
        ↓
frozen-composition characteristic slope limiting
        ↓
MUSCL-Hancock face prediction
        ↓
general-EOS Rusanov or HLLC flux
        ↓
conservative finite-volume update
```

Species face fluxes are corrected so their sum equals the total mass flux to roundoff.

## Reaction-flow coupling

The first coupled integrator uses Strang splitting:

```text
reaction(dt/2)
      ↓
hydro(dt)
      ↓
reaction(dt/2)
```

Each cell reaction solve holds density, all momentum components, and total-energy density fixed. Composition changes are written back to `rho*Y_k`, and temperature is recovered from the unchanged specific internal energy. This avoids adding a separate heat-release source on top of formation-energy-inclusive NASA7 internal energy.

## Verification separation

The architecture retains independent gates for:

1. constant-`gamma` hydro;
2. passive multispecies transport;
3. NASA7 thermodynamics;
4. zero-dimensional elementary chemistry;
5. composition-dependent reactive hydro;
6. reaction-flow splitting.

The homogeneous reactive field must reduce to independent zero-dimensional cell chemistry. The nonuniform hotspot must create finite pressure and velocity responses while preserving global mass, momentum, and total energy.

## Design constraints

1. Lower-order and constant-`gamma` baselines remain active in CI.
2. Invalid density, pressure, composition, temperature, or energy inversion fails explicitly.
3. Formation-energy offsets are retained in `h` and `u`.
4. Reverse rates use the same NASA7 records as the energy equation.
5. Species fluxes close exactly to the shared mass flux.
6. Chemistry does not independently modify `rhoE` in the adiabatic constant-volume substep.
7. The current characteristic basis assumes frozen composition across each acoustic solve.
8. Rusanov remains the robustness baseline. HLLC is the verified
   contact-resolving general-EOS intermediate; it is not labeled as PeleC
   Riemann parity.
9. The current four-reaction chemistry subset is not a complete H2/O2 mechanism.

## Reactive PPM path

`reactive_1d_mod` now keeps PCM, characteristic PLM, and monotone primitive PPM as separate reconstruction paths. PPM uses a wider explicitly filled stencil and SSPRK3, while the Riemann solver remains selected independently.
