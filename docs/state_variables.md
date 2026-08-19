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

with a separate scalar temperature and density. The fast H2/O2 subset fixes the ordering

```text
1 H2
2 H
3 O
4 O2
5 OH
6 H2O
7 N2
```

The full pressure-dependent mechanism uses

```text
1 H2
2 H
3 O
4 O2
5 OH
6 H2O
7 HO2
8 H2O2
9 AR
10 N2
```

and associates each entry with a NASA7 record and molecular weight. Concentrations and production rates use

```text
C_k        kmol/m^3
wdot_k     kmol/(m^3 s)
dY_k/dt    1/s
```

The reaction record stores two runtime stoichiometric vectors, a reaction-kind tag, high- and low-pressure Arrhenius records, a third-body efficiency vector, optional Troe parameters, and a reversible flag. Reactor density is fixed by the initial `(T,p,Y)` state; temperature is derived from the fixed specific internal energy rather than advanced as an independent ODE component.

This separation is intentional: `0.10.0` chemistry evidence does not imply that current hydro pressure, sound speed, or temperature are composition dependent.
