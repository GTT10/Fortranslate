# PeleF architecture

## Executable split

PeleF exposes four serial verification drivers over shared numerical modules.

```text
pelef
  └─ one-dimensional PCM / PLM / PeleC-style Godunov paths

pelef2d
  └─ two-dimensional periodic CTU-style Euler path

pelef_ms
  └─ passive runtime-multispecies transport over the 1D/2D hydro core

pelef0d
  └─ NASA7 mixture thermodynamics and a constant-volume toy reactor
```

The split prevents an unverified general-EOS or chemistry path from silently changing the already pinned constant-`gamma` hydro results.

## Hydro state and physics boundary

The verified hydro prefix is

```text
(rho, rho*u, rho*v, rho*w, rho*E).
```

The multispecies path appends conserved `rho*Y_k` components. Its hydrodynamic pressure and sound speed still come from the original constant-`gamma` EOS. The new NASA7 modules are therefore a separate thermodynamic layer rather than an implicit replacement inside the Riemann solvers.

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
        ↓
isomerization_reactor_mod
  └─ first-order A -> B source integration at constant volume
```

`temperature_from_internal_energy` uses the common NASA7 validity interval of all participating species. Newton updates use mixture `cv`; any step leaving the current bracket is replaced by bisection.

## Zero-dimensional reactor path

```text
reactor namelist
  ↓
synthetic equal-molecular-weight NASA7 species A and B
  ↓
initial specific internal energy
  ↓
RK4 composition stages
  ├─ isothermal: fixed temperature
  └─ adiabatic: solve e(Y,T) = e_initial at every stage
  ↓
CSV output and independent monotonicity/energy checks
```

The synthetic species isolate thermodynamic and integration algebra. They do not represent a physical fuel mechanism.

## One- and two-dimensional hydro paths

The existing hydro architecture remains unchanged:

- 1D: PCM, componentwise PLM, or characteristic PLM; Rusanov or qualified PeleC Riemann; SSPRK2 or time-centered Godunov update.
- 2D: directional rotation, provisional fluxes, transverse half-step corrections, final Riemann fluxes, and one unsplit conservative update.
- Multispecies: species fluxes consume the same face mass flux as hydro and satisfy `sum_k F_(rho Y_k) = F_rho`.

## Design constraints

1. Constant-`gamma` hydro baselines remain bitwise isolated from the new thermodynamics layer.
2. NASA7 records carry explicit validity bounds and molecular weights.
3. Invalid composition, out-of-range temperature, or unbracketed internal energy fails explicitly.
4. Formation-energy offsets are retained in `h` and `u`; they are not discarded as merely sensible heat.
5. Adiabatic reactor stages recover temperature from the same fixed internal-energy target.
6. A chemistry subsystem is not described as PeleC-compatible until a real mechanism, reaction rates, and an external parity case are present.
7. General-EOS hydro coupling must receive its own Riemann, CFL, and conservation regressions rather than reusing the constant-`gamma` evidence by assumption.
