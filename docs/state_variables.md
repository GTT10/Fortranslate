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

## Planned multispecies extension

Phase 3 will append conserved species densities `rho*Y_k` after the current Euler state. The five-variable layout must remain the zero-species specialization so existing hydro kernels and regressions continue to compile and run unchanged.
