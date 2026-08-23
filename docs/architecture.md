# PeleF architecture

## Executable split

PeleF exposes nine serial verification drivers and five optional MPI drivers
over shared numerical and physical-property modules.

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

pelef0d_h2o2_full
  └─ full pressure-dependent H2/O2 kinetics and an implicit reactor

pelef_transport_probe
  └─ qualified dilute-gas mixture transport coefficients

pelef_reactive_1d
  └─ NASA7 general-EOS reactive Euler with PLM/PPM, HLLC, and Strang splitting

pelef_reactive_2d
  └─ NASA7 reactive Euler, physical boundaries, transport, and CTU correction
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


## Shared molecular transport

```text
transport_database_mod
  └─ pinned Lennard-Jones records for the seven-species subset
        ↓
mixture_transport_mod
  ├─ Chapman--Enskog pure viscosity
  ├─ Wilke mixture viscosity
  ├─ modified-Eucken / Mathur thermal conductivity
  ├─ Chapman--Enskog binary diffusion
  └─ mixture-averaged species diffusion
        ↓
reactive_diffusive_flux_x
  ├─ Newtonian viscous stress and viscous work
  ├─ Fourier heat flux
  ├─ mole-fraction and optional pressure diffusion driving forces
  ├─ correction velocity enforcing sum(j_k) = 0
  └─ species-enthalpy diffusion energy flux
```

Transport data are pinned to the same Cantera `h2o2.yaml` provenance as the
current thermodynamics/chemistry subset. The implemented coefficient model is
a dilute ideal-gas subset, not the full PelePhysics generated polynomial
transport layer.

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
reconstruction selector
  ├─ PCM
  ├─ characteristic PLM + MUSCL-Hancock tracing
  ├─ componentwise monotone PPM + SSPRK3
  └─ characteristic PPM profile integration
       ├─ optional contact steepening
       └─ optional shock flattening
        ↓
general-EOS Rusanov or HLLC flux
        ↓
conservative finite-volume update

optional transport branch
  ├─ face-centered viscous / conductive / species fluxes
  ├─ explicit SSPRK2 diffusion update
  └─ parabolic dx^2 / diffusivity timestep gate
```

Advective species face fluxes close to the total mass flux. Diffusive species
fluxes use a correction velocity so their sum is zero to roundoff.

## Reactive two-dimensional CTU path

The reactive 2D state is stored as `state(variable,nx,ny)` with a synchronized
`temperature(nx,ny)` field. The normal predictor is selected independently from
the Riemann solver and supports PCM, frozen-composition characteristic PLM, or
time-traced characteristic PPM. For the y direction, momentum and primitive
velocity components are rotated into the x-normal ordering, evaluated by the
same predictor and HLLC/Rusanov kernels, and rotated back.

```text
cell-centered conserved state
        ↓
NASA7 conserved-to-primitive recovery
        ↓
normal predictor selector
  ├─ PCM
  ├─ characteristic PLM
  └─ characteristic PPM
       ├─ five-point parabolic reconstruction
       ├─ u-c / u / u+c profile integration
       ├─ optional bounded contact steepening
       └─ optional shock flattening
        ↓
provisional x/y HLLC fluxes
        ↓
conservative transverse half-step correction
  U_face* = U_face - dt/(2 d_t) (F_t,hi - F_t,lo)
        ↓
EOS/positivity bisection on the complete conserved face state
        ↓
final directional HLLC fluxes
        ↓
unsplit two-dimensional conservative update
```

The transverse limiter acts on the complete conserved vector, including every
species density and total energy. Because each directional species-flux block
closes to the corresponding mass flux, the corrected face state retains
`sum(rho*Y_k)=rho` rather than repairing species independently after the hydro
correction. Both x-normal and y-normal characteristic-PLM/PPM reductions agree
with the corresponding 1D update at roundoff.

Chemistry uses the same cell-local adiabatic constant-volume solver as the 1D
path and is Strang split around the unsplit CTU hydro step. This path currently
requires periodic boundaries. The characteristic-PPM option supplies a
PeleC-style normal predictor in each coordinate direction before the existing
full-state CTU correction. Complete PeleC multidimensional PPM transverse/corner
tracing remains intentionally outside the current claim.

## Reaction-flow coupling

With molecular transport disabled, the coupled integrator retains the
reaction--hydro--reaction Strang sequence. With transport enabled, the symmetric
composition is:

```text
reaction(dt/2)
      ↓
transport(dt/2)
      ↓
hydro(dt)
      ↓
transport(dt/2)
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
6. reaction-flow splitting;
7. reactive two-dimensional CTU and dimensional reduction;
8. one- and two-dimensional molecular transport and Cantera coefficient
   qualification;
9. physical boundaries and full pressure-dependent H2/O2 chemistry;
10. MPI decomposition, multispecies hydro, implicit chemistry, molecular
    transport, and coupled reactive splitting.

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
9. The four-reaction chemistry subset remains a lightweight regression path;
   the selectable ten-species, 29-reaction mechanism is the full H2/O2 path.
10. Contact steepening is explicitly bounded to half of the canonical detector
    strength until a complete general-EOS PPM/HLLC characteristic system is
    available.
11. The present transport layer excludes Soret, Dufour, multicomponent Stefan--
    Maxwell diffusion, polar corrections, and bulk viscosity.
12. Molecular transport is qualified in serial 1D/2D and distributed 1D paths;
    Soret and multicomponent diffusion remain excluded.

## Reactive PPM path

`reactive_1d_mod` keeps four independently selectable paths:

- `pcm`, the first-order robustness baseline;
- `characteristic_plm`, the frozen-composition MUSCL-Hancock path;
- `ppm`, the semidiscrete componentwise monotone path advanced by SSPRK3;
- `characteristic_ppm`, a time-centered normal predictor using PeleC's
  five-point parabolic reconstruction and `u-c`, `u`, `u+c` profile
  integration.

The characteristic PPM path carries species and transverse velocities on the
middle wave, projects density/normal velocity/pressure over the frozen mixture
acoustic basis, and converts the final face state through the NASA7 EOS. Its
one-dimensional shock-flattening coefficient follows PeleC `Godunov.H`.
Contact steepening is a separate Colella--Woodward-style detector applied to
density and species only. Both controls are opt-in and are rejected by the
configuration reader for other reconstruction modes.


## Reactive two-dimensional molecular transport

`reactive_transport_2d_mod` evaluates x/y face transport fluxes and advances
their conservative divergence independently of the CTU hydro operator.


## Physical boundary layer

`reactive_boundary_2d_mod` owns four typed faces and samples periodic, wall,
inflow, or outflow ghost states. Hydro and molecular transport use explicit
lower/upper face arrays. Solid walls receive a pressure-only inviscid flux,
mirrored velocity/temperature transport gradients, and zero species flux.

## Full pressure-dependent H2/O2 chemistry

The reactive applications dispatch either the seven-species elementary model or a ten-species, 29-reaction model. The full path reuses the variable-width conserved state and advances each cell with the implicit constant-volume reactor.

## MPI one-dimensional path

`mpi_domain_1d_mod` owns uneven contiguous decomposition, periodic nonblocking
halo exchange, global reductions, and ordered `MPI_Gatherv` output. The MPI
drivers allocate only rank-local state plus two ghost cells.

```text
rank-local conserved state + temperature
        ↓
periodic state/temperature halo exchange
        ↓
global hydro/transport timestep reduction
        ↓
implicit chemistry(dt/2)
        ↓
SSPRK2 molecular transport(dt/2)
        ↓
conservative general-EOS Rusanov hydro(dt)
        ↓
SSPRK2 molecular transport(dt/2)
        ↓
implicit chemistry(dt/2)
```

Every coupled attempt is transactional. A failure on any rank is reduced across
the communicator, all ranks restore the pre-attempt state, and the scheduler
retries a smaller interval. CI compares complete gathered fields for 1, 2, 4,
and 8 ranks in both Debug and Release builds.
