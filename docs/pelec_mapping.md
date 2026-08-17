# PeleC to PeleF mapping

This table records responsibility mapping, not line-by-line source translation.

| PeleC reference | PeleF implementation | Current status |
|---|---|---|
| `Source/main.cpp` | `app/pelef.F90` | Minimal serial driver implemented |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base conserved/primitive indices implemented |
| PelePhysics ideal-gas EOS calls | `src/physics/eos_ideal_mod.F90` | Constant-gamma ideal gas implemented |
| Conserved/primitive conversion in hydro kernels | `src/hydro/state_conversion_mod.F90` | Implemented |
| `Source/Riemann.H` local Lax–Friedrichs path | `src/hydro/riemann_rusanov_mod.F90` | Independent Rusanov implementation completed |
| `Source/Godunov.*` flux-divergence responsibility | `src/hydro/finite_volume_mod.F90` | First-order 1D operator implemented |
| `Source/Timestep.H` hydro CFL responsibility | `src/driver/time_integrator_mod.F90` | Implemented for 1D ideal gas |
| `Source/Advance.cpp` time-advance responsibility | `src/driver/time_integrator_mod.F90` | SSPRK2 implemented; PeleC MOL/SDC parity not yet attempted |
| `Exec/RegTests/Sod` | `cases/sod`, `tools/compare_sod.py` | Exact-solution regression implemented |
| `Source/PLM.H` | future `src/hydro/reconstruction_plm_mod.F90` | Not started |
| `Source/PPM.*` | future `src/hydro/reconstruction_ppm_mod.F90` | Not started |
| `Source/WENO.H` | future `src/hydro/reconstruction_weno_mod.F90` | Not started |
| `Source/Diffusion.*` | future `src/diffusion/` | Not started |
| `Source/React.cpp` | future `src/chemistry/` | Not started |
| `Source/PeleCAmr.*` | future `src/amr/` | Not started |
| `Source/EB.*` | future `src/eb/` | Not started |
| `Source/LES.*` | future `src/les/` | Not started |
| `Source/Particle.cpp` | future `src/particles/` | Not started |

## Traceability rule

A row moves to “implemented” only when the Fortran subsystem has an automated test. A name or file existing without numerical validation does not count as a completed port.
