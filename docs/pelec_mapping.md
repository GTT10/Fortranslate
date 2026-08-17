# PeleC to PeleF mapping

This table maps responsibilities, not lines of C++ to lines of Fortran.

| PeleC reference | PeleF implementation | Status |
|---|---|---|
| `Source/main.cpp` | `app/pelef.F90` | Minimal serial driver |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base state indices |
| PelePhysics ideal-gas calls | `src/physics/eos_ideal_mod.F90` | Constant-`gamma` subset |
| state conversion kernels | `src/hydro/state_conversion_mod.F90` | Implemented and tested |
| AMReX boundary fills | `src/hydro/boundary_conditions_mod.F90` | 1D outflow/periodic subset |
| local Lax-Friedrichs path | `src/hydro/riemann_rusanov_mod.F90` | Implemented |
| `Source/Riemann.H` | `src/hydro/riemann_pelec_mod.F90` | Qualified single-species ideal-gas subset |
| Riemann selection | `src/hydro/riemann_flux_mod.F90` | Explicit Rusanov/PeleC dispatch |
| `Source/PLM.H` slope scaffold | `slope_limiter_mod.F90`, `reconstruction_plm_mod.F90` | Componentwise baseline |
| `Source/PLM.H` characteristic projection/tracing | `reconstruction_pelec_plm_mod.F90` | Qualified 1D constant-`gamma` subset |
| `Source/Godunov.*` update responsibility | `finite_volume_mod.F90`, `time_integrator_mod.F90` | PCM/PLM MOL and traced Godunov paths |
| `Source/Timestep.H` | `time_integrator_mod.F90` | 1D hydro CFL subset |
| `Exec/RegTests/Sod` | `cases/sod`, `tools/compare_sod.py` | Multiple exact-solution gates |
| `Exec/RegTests/Shu-Osher` | `cases/shu_osher`, `tools/check_shu_osher.py` | Deterministic signature gate |
| `Source/PPM.*` | future `src/hydro/reconstruction_ppm_mod.F90` | Not started |
| `Source/WENO.H` | future `src/hydro/reconstruction_weno_mod.F90` | Not started |
| `Source/Diffusion.*` | future `src/diffusion/` | Not started |
| `Source/React.cpp` | future `src/chemistry/` | Not started |
| `Source/PeleCAmr.*` | future `src/amr/` | Not started |
| `Source/EB.*` | future `src/eb/` | Not started |
| `Source/LES.*` | future `src/les/` | Not started |
| `Source/Particle.cpp` | future `src/particles/` | Not started |

A row is marked implemented only when its Fortran subsystem has an automated numerical gate. The `Riemann.H` and `PLM.H` rows remain explicitly qualified because general EOS, species, flattening, EB, and multidimensional logic are absent.
