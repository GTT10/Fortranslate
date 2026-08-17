# State variables

## Phase-1 conserved array

The solver advances five Euler conserved variables using Fortran-native one-based indices:

| PeleF index | Symbol | Meaning |
|---:|---|---|
| `IRHO = 1` | `rho` | mass density |
| `IMX = 2` | `rho*u` | x-momentum density |
| `IMY = 3` | `rho*v` | y-momentum density |
| `IMZ = 4` | `rho*w` | z-momentum density |
| `IET = 5` | `rho*E` | total-energy density |

Although the first solver is one-dimensional, transverse momentum components are retained so the state layout extends naturally to multidimensional fluxes.

## Reserved base-state positions

`IEI = 6` and `ITEM = 7` reserve the planned PeleF base-state positions for internal-energy density and temperature. They are derived quantities in Phase 1 and are not independently advanced. This avoids storing thermodynamically inconsistent duplicate values in the minimal solver.

## Primitive array

| PeleF index | Meaning |
|---:|---|
| `QRHO = 1` | density |
| `QU = 2` | x velocity |
| `QV = 3` | y velocity |
| `QW = 4` | z velocity |
| `QP = 5` | pressure |
