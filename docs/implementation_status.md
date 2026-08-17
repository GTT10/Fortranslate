# Implementation status

## Phase 0 — project and test infrastructure

- [x] CMake build
- [x] GNU Fortran support
- [x] Debug bounds/FPE checks
- [x] CTest harness
- [x] GitHub Actions workflow
- [x] architecture, mapping, state, and parity documents

## Phase 1 — minimal compressible Euler solver

- [x] uniform one-dimensional mesh
- [x] constant-gamma ideal-gas EOS
- [x] primitive/conserved conversion
- [x] outflow ghost-cell boundaries
- [x] Rusanov flux
- [x] conservative finite-volume divergence
- [x] hydro CFL timestep
- [x] SSPRK2 time integration
- [x] Sod input and initialization
- [x] CSV output
- [x] exact Sod comparison
- [ ] Shu–Osher regression
- [ ] Sedov regression
- [ ] smooth-advection regression
- [ ] isentropic-vortex regression

## Next implementation slice

Add piecewise-linear reconstruction with a selectable limiter while retaining the existing first-order path as a regression baseline. Validate convergence and shock behavior before beginning the PeleC-style approximate Riemann solver.
