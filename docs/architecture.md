# PeleF architecture

## Executable split

PeleF exposes five serial verification drivers over shared numerical and physical-property modules.

```text
pelef
  └─ one-dimensional PCM / PLM / PeleC-style Godunov paths

pelef2d
  └─ two-dimensional periodic CTU-style Euler path

pelef_ms
  └─ passive runtime-multispecies transport over the 1D/2D hydro core

pelef0d
  └─ NASA7 mixture thermodynamics and a synthetic isomerization reactor

pelef0d_h2o2
  └─ generated reversible elementary kinetics and an H2/O2 reactor
```

The split prevents an unverified general-EOS or chemistry path from silently changing the already pinned constant-`gamma` hydro results.

## Hydro state and physics boundary

The verified hydro prefix is

```text
(rho, rho*u, rho*v, rho*w, rho*E).
```

The multispecies path appends conserved `rho*Y_k` components. Its hydrodynamic pressure and sound speed still come from the original constant-`gamma` EOS. NASA7 and chemistry modules remain an explicit, independent layer until a later general-EOS hydro milestone adds new Riemann, characteristic, CFL, and conservation evidence.

## Thermodynamics layer

```text
thermo_database_mod
  └─ named NASA7 species records and molecular weights
        ↓
nasa7_thermo_mod
  └─ species cp, cv, h, u, and standard-state entropy
        ↓
mixture_thermo_mod
  ├─ mass-based mixture molecular weight and gas constant
  ├─ mixture cp, cv, gamma, h, and u
  ├─ ideal-gas pressure, density, and frozen sound speed
  └─ safeguarded e -> T inversion
```

`temperature_from_internal_energy` uses the common NASA7 validity interval of all participating species. Newton updates use mixture `cv`; any step leaving the current bracket is replaced by bisection.

## Generated elementary-kinetics path

```text
mechanisms/h2o2_elementary.json
        ↓
tools/generate_elementary_mechanism.py
        ↓
src/generated/h2o2_elementary_mechanism_mod.F90
        ↓
elementary_kinetics_mod
  ├─ arbitrary stoichiometric vectors
  ├─ forward Arrhenius rate
  ├─ NASA7 equilibrium constant and reverse rate
  ├─ progress rates
  └─ species molar production rates
        ↓
constant_volume_reactor_mod
  ├─ dY/dt from W_k * wdot_k / rho
  ├─ adaptive RK4 step doubling
  └─ stage-wise e(Y,T) = e_initial inversion
        ↓
pelef0d_h2o2
  └─ CSV history, structural checker, and Cantera parity
```

The normalized JSON file is the authoring format for the current generated subset. The generated Fortran source is committed, and CI regenerates it into a temporary file and requires byte-for-byte agreement. The companion Cantera YAML expresses the same species and four reactions for an independent runtime reference.

## Zero-dimensional reactor paths

The synthetic reactor remains a small algebraic gate:

```text
A -> B
```

The H2/O2 reactor is the first physical chemistry gate. It uses seven species and four reversible elementary reactions. Both reactors are adiabatic and constant volume when energy coupling is enabled; temperature is recovered from the same fixed initial specific internal energy at every accepted or trial stage.

## One- and two-dimensional hydro paths

The existing hydro architecture remains unchanged:

- 1D: PCM, componentwise PLM, or characteristic PLM; Rusanov or qualified PeleC Riemann; SSPRK2 or time-centered Godunov update.
- 2D: directional rotation, provisional fluxes, transverse half-step corrections, final Riemann fluxes, and one unsplit conservative update.
- Multispecies: species fluxes consume the same face mass flux as hydro and satisfy `sum_k F_(rho Y_k) = F_rho`.

## Design constraints

1. Constant-`gamma` hydro baselines remain isolated from NASA7 thermodynamics and chemistry.
2. NASA7 records carry explicit validity bounds and molecular weights.
3. Invalid composition, temperature, stoichiometry, or unbracketed internal energy fails explicitly.
4. Formation-energy offsets are retained in `h` and `u`.
5. Reverse rates are derived from the same NASA7 records used by the reactor energy equation.
6. Generated mechanism source must be reproducible from the committed normalized mechanism input.
7. Production-rate parity is evaluated at the exact PeleF state, separately from trajectory parity, so kinetic-kernel errors are not confused with integrator divergence.
8. Explicit RK4 evidence is limited to this non-stiff verification case; it is not a substitute for a stiff chemistry integrator.
9. General-EOS hydro coupling must receive its own Riemann, CFL, and conservation regressions.
