# PeleC to PeleF mapping

This table maps responsibilities, not source lines.

| PeleC reference | PeleF implementation | Status |
|---|---|---|
| `Source/main.cpp` | `app/pelef.F90` | Minimal serial problem driver |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base single-species state indices |
| PelePhysics ideal-gas calls | `src/physics/eos_ideal_mod.F90` | Constant-`gamma` ideal gas |
| AMReX boundary fill | `src/hydro/boundary_conditions_mod.F90` | 1D outflow and periodic fills |
| `Source/Riemann.H` LF path | `src/hydro/riemann_rusanov_mod.F90` | Independent Rusanov implementation |
| `Source/Riemann.H` acoustic solver | `src/hydro/riemann_pelec_mod.F90` | Qualified single-species ideal-gas subset |
| `Source/PLM.H` characteristic projection/tracing | `src/hydro/reconstruction_pelec_plm_mod.F90` | Qualified 1D subset |
| `Source/PLM.H::plm_slope` | `pelec_limited_slope` | Order 2 and 4 regular-cell formulas verified |
| `Source/Godunov.H::flatten` | `pelec_flattening_coefficient` | 1D regular-cell formula verified |
| `Source/Godunov.*` conservative update | `src/hydro/finite_volume_mod.F90`, `src/driver/time_integrator_mod.F90` | 1D time-centered update |
| `Exec/RegTests/Sod` | `cases/sod`, `tools/compare_sod.py` | Exact-solution regressions |
| `Exec/RegTests/Shu-Osher` | `cases/shu_osher`, `tools/check_shu_osher.py` | Deterministic shock-wave interaction gate |
| `Exec/RegTests/Sedov` | `cases/sedov`, `tools/check_sedov.py` | Independent planar strong-blast gate |
| `Source/PPM.*` | future `src/hydro/reconstruction_ppm_mod.F90` | Not started |
| multidimensional Godunov corrections | future 2D hydro modules | Not started |
| `Source/Diffusion.*` | future `src/diffusion/` | Not started |
| `Source/React.cpp` | future `src/chemistry/` | Not started |
| `Source/PeleCAmr.*` | future `src/amr/` | Not started |
| `Source/EB.*` | future `src/eb/` | Not started |
| `Source/LES.*` | future `src/les/` | Not started |
| `Source/Particle.cpp` | future `src/particles/` | Not started |

A row is called implemented only when its Fortran subsystem has an automated numerical gate. The current PLM, flattening, and Riemann rows remain explicitly qualified because general EOS, species, embedded boundaries, and multidimensional behavior are outside their tested scope.
