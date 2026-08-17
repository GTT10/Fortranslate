# PeleC to PeleF mapping

This table records responsibility mapping, not line-by-line source translation.

| PeleC reference | PeleF implementation | Current status |
|---|---|---|
| `Source/main.cpp` | `app/pelef.F90` | Serial problem-selecting driver implemented |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base conserved/primitive indices implemented |
| PelePhysics ideal-gas EOS calls | `src/physics/eos_ideal_mod.F90` | Constant-gamma ideal gas implemented |
| Conserved/primitive conversion in hydro kernels | `src/hydro/state_conversion_mod.F90` | Implemented |
| AMReX boundary-fill responsibility | `src/hydro/boundary_conditions_mod.F90` | 1D outflow and periodic fills implemented |
| `Source/Riemann.H` local Lax-Friedrichs path | `src/hydro/riemann_rusanov_mod.F90` | Independent Rusanov implementation completed |
| `Source/Riemann.H` acoustic star-state path | `src/hydro/riemann_pelec_mod.F90` | Single-species constant-gamma reduction verified; general EOS/species parity pending |
| Riemann-method selection | `src/hydro/riemann_flux_mod.F90` | Explicit Rusanov/PeleC-style dispatch verified |
| `Source/Godunov.*` flux-divergence responsibility | `src/hydro/finite_volume_mod.F90` | Selectable first/second-order 1D operator implemented |
| `Source/PLM.H` | `src/hydro/slope_limiter_mod.F90`, `src/hydro/reconstruction_plm_mod.F90` | Componentwise primitive PLM scaffold verified; characteristic parity pending |
| `Source/Timestep.H` hydro CFL responsibility | `src/driver/time_integrator_mod.F90` | Implemented for 1D ideal gas |
| `Source/Advance.cpp` time-advance responsibility | `src/driver/time_integrator_mod.F90` | SSPRK2 implemented; PeleC MOL/SDC parity not yet attempted |
| `Exec/RegTests/Sod` | `cases/sod`, `tools/compare_sod.py` | PCM, PLM/Rusanov, and PLM/PeleC-style exact-solution regressions implemented |
| `Exec/RegTests/Shu-Osher` | `cases/shu_osher`, `tools/check_shu_osher.py` | Deterministic shock-density-wave regression implemented |
| Smooth-method verification | `tests/regression/test_entropy_wave*.F90` | Second-order gates for both flux choices implemented |
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

Qualified rows state exactly which subset is covered. The current Riemann implementation reproduces the ideal-gas acoustic star-state and wave interpolation structure, but not PelePhysics general-EOS/species behavior, rotating-frame energy, or multidimensional Godunov coupling.
