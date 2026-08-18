# State variables

## Conserved state

Both the one-dimensional and two-dimensional solvers advance five Euler conserved variables using Fortran-native one-based indices:

| PeleF index | Symbol | Meaning |
|---:|---|---|
| `IRHO = 1` | `rho` | mass density |
| `IMX = 2` | `rho*u` | x-momentum density |
| `IMY = 3` | `rho*v` | y-momentum density |
| `IMZ = 4` | `rho*w` | z-momentum density |
| `IET = 5` | `rho*E` | total-energy density |

The transverse momentum components were retained from the first one-dimensional milestone and are now used directly by the 2D directional flux and CTU operators.

## Array layouts

The verified 1D drivers use

```fortran
state(variable, 0:nx+1)
```

with one conserved ghost cell on each side. Wider reconstruction stencils are built in temporary primitive arrays.

The periodic 2D driver uses

```fortran
state(variable, 1:nx, 1:ny)
```

for interior cells only. Periodic neighbors are obtained with wrapped indices, and one shared flux is stored for each right and upper cell face.

## Reserved base-state positions

`IEI = 6` and `ITEM = 7` reserve planned base-state positions for internal-energy density and temperature. They are derived quantities and are not independently advanced in the constant-`gamma` solver.

## Primitive state

| PeleF index | Meaning |
|---:|---|
| `QRHO = 1` | density |
| `QU = 2` | x velocity |
| `QV = 3` | y velocity |
| `QW = 4` | z velocity |
| `QP = 5` | pressure |

For a y-normal Riemann problem, `QU` and `QV` are exchanged before calling the x-direction solver and exchanged again on return.

## Runtime multispecies layout

The passive multispecies path uses

```text
1:5              rho, rho*u, rho*v, rho*w, rho*E
6                derived rho*e
7                derived temperature proxy
8:7+nspecies     rho*Y_1 ... rho*Y_N
```

`multispecies_nvar(nspecies) = 7 + nspecies`, and `species_component(k) = 7 + k`. Positions 6 and 7 are synchronized from the conserved hydro prefix after every update and are never evolved by an independent flux. The five-variable hydro layout remains unchanged for existing single-species solvers.

Every accepted multispecies state must satisfy non-negative species densities and

```text
sum_k rho*Y_k = rho
```

within a scaled tolerance. Invalid closure is rejected rather than silently renormalizing the cell state. Face mass fractions are bounded and normalized before constructing species fluxes.


## Thermodynamics-only composition state

The NASA7 and zero-dimensional reactor modules do not reuse the hydro-derived temperature proxy. They take an explicit mass-fraction vector `Y(:)` and temperature, or recover temperature from a target specific internal energy. Species molecular weights and polynomial records live outside the conserved flow array.

This separation is intentional: the `0.8.0` thermodynamics evidence does not imply that the current constant-`gamma` hydro fluxes are composition dependent. A later general-EOS hydro milestone must replace the derived pressure, sound-speed, and temperature path explicitly and add new conservation/parity gates.
