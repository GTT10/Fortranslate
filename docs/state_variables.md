# State variables

## Conserved hydro state

Both the one-dimensional and two-dimensional solvers advance five Euler conserved variables using Fortran-native one-based indices:

| PeleF index | Symbol | Meaning |
|---:|---|---|
| `IRHO = 1` | `rho` | mass density |
| `IMX = 2` | `rho*u` | x-momentum density |
| `IMY = 3` | `rho*v` | y-momentum density |
| `IMZ = 4` | `rho*w` | z-momentum density |
| `IET = 5` | `rho*E` | total-energy density |

The verified 1D drivers use `state(variable,0:nx+1)`. The periodic 2D driver uses `state(variable,1:nx,1:ny)` and wrapped neighbor indices.

## Primitive state

| PeleF index | Meaning |
|---:|---|
| `QRHO = 1` | density |
| `QU = 2` | x velocity |
| `QV = 3` | y velocity |
| `QW = 4` | z velocity |
| `QP = 5` | pressure |

For a y-normal Riemann problem, `QU` and `QV` are exchanged before calling the x-direction solver and exchanged again on return.

## Runtime passive-multispecies layout

The passive multispecies path uses

```text
1:5              rho, rho*u, rho*v, rho*w, rho*E
6                derived rho*e
7                derived temperature proxy
8:7+nspecies     rho*Y_1 ... rho*Y_N
```

Positions 6 and 7 are synchronized from the conserved hydro prefix and are not independently evolved. Every accepted state must satisfy non-negative species densities and

```text
sum_k rho*Y_k = rho
```

within a scaled tolerance.

## Thermodynamics and chemistry composition state

NASA7 and reactor modules use an explicit mass-fraction vector

```fortran
Y(nspecies)
```

with a separate scalar temperature and density. The H2/O2 generated subset fixes the ordering

```text
1 H2
2 H
3 O
4 O2
5 OH
6 H2O
7 N2
```

and associates each entry with a NASA7 record and molecular weight. Concentrations and production rates use

```text
C_k        kmol/m^3
wdot_k     kmol/(m^3 s)
dY_k/dt    1/s
```

The reaction record stores two runtime stoichiometric vectors of length `nspecies` plus an Arrhenius triplet and reversible flag. Reactor density is fixed by the initial `(T,p,Y)` state; temperature is derived from the fixed specific internal energy rather than advanced as an independent ODE component.

This explicit-vector representation remains the zero-dimensional chemistry interface. The older `pelef_ms` path is still passive and constant-`gamma`; the separate `pelef_reactive_1d` layout below couples the same NASA7 composition to pressure, temperature, sound speed, and reaction-flow splitting.

## Reactive one-dimensional conserved layout

The `pelef_reactive_1d` path uses a compact runtime layout without the passive path's derived slots:

```text
1:5              rho, rho*u, rho*v, rho*w, rho*E
6:5+nspecies     rho*Y_1 ... rho*Y_N
```

For the current seven-species mechanism,

```text
6  rho*Y_H2
7  rho*Y_H
8  rho*Y_O
9  rho*Y_O2
10 rho*Y_OH
11 rho*Y_H2O
12 rho*Y_N2
```

The corresponding primitive vector is

```text
rho, u, v, w, p, Y_H2, Y_H, Y_O, Y_O2, Y_OH, Y_H2O, Y_N2.
```

Temperature and frozen sound speed are derived from `(rho, rhoE, momentum, rhoY)` by removing kinetic energy, recovering mass fractions, solving `u(Y,T)=e`, and evaluating the NASA7 mixture EOS. They are cached only as synchronized guesses/diagnostics and are never independently conserved.

During a chemistry substep, `rho`, momentum, and `rhoE` remain fixed while `rhoY_k` changes. During the hydro substep, every conserved component is advanced by a face-flux divergence. Every accepted cell must satisfy

```text
rho > 0
p > 0
T inside the common NASA7 range
Y_k >= 0
sum_k Y_k = 1.
```
