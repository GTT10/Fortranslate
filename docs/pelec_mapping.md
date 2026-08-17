# PeleC to PeleF mapping

This table maps responsibilities, not source lines.

| PeleC reference | PeleF implementation | Status |
|---|---|---|
| `Source/main.cpp` | `app/pelef.F90`, `app/pelef2d.F90` | Separate serial 1D and 2D drivers |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base single-species state indices |
| PelePhysics ideal-gas calls | `src/physics/eos_ideal_mod.F90` | Constant-`gamma` ideal gas |
| AMReX mesh/geometry responsibility | `mesh_mod`, `mesh_2d_mod` | Uniform 1D and Cartesian 2D meshes |
| AMReX boundary fill | `boundary_conditions_mod` and periodic index wrapping in `ctu_2d_mod` | 1D outflow/periodic and 2D periodic subset |
| `Source/Riemann.H` LF path | `riemann_rusanov_mod` | Independent Rusanov implementation |
| `Source/Riemann.H` acoustic solver | `riemann_pelec_mod` | Qualified single-species ideal-gas subset |
| Direction-dependent flux assembly | `directional_flux_mod` | x/y momentum rotation and y fluxes verified |
| `Source/PLM.H` characteristic projection/tracing | `reconstruction_pelec_plm_mod` | Qualified 1D regular-cell subset |
| `Source/PLM.H::plm_slope` | `pelec_limited_slope` | Order 2 and 4 formulas verified |
| `Source/Godunov.H::flatten` | `pelec_flattening_coefficient` | 1D regular-cell formula verified |
| `Source/Godunov.*` transverse update responsibility | `src/hydro/ctu_2d_mod.F90` | Qualified periodic 2D regular-grid CTU-style subset |
| 1D conservative update | `finite_volume_mod`, `time_integrator_mod` | SSPRK2 and time-centered Godunov paths |
| 2D conservative update | `ctu_2d_mod` | Normal prediction, provisional fluxes, transverse correction, final update |
| `Exec/RegTests/Sod` | `cases/sod`, `tools/compare_sod.py` | Exact-solution regressions |
| `Exec/RegTests/Shu-Osher` | `cases/shu_osher`, `tools/check_shu_osher.py` | Deterministic shock-wave gate |
| `Exec/RegTests/Sedov` | `cases/sedov`, `tools/check_sedov.py` | Independent planar strong-blast gate |
| multidimensional smooth verification | `cases/isentropic_vortex`, `test_isentropic_vortex_2d`, `check_isentropic_vortex.py` | Analytical 2D convergence and app gate |
| `Exec/RegTests/MultiSpecSod` | future multispecies state and case | Not started |
| `Source/PPM.*` | future `reconstruction_ppm_mod` | Not started |
| `Source/WENO.H` | future `reconstruction_weno_mod` | Not started |
| `Source/Diffusion.*` | future `src/diffusion/` | Not started |
| `Source/React.cpp` | future `src/chemistry/` | Not started |
| `Source/PeleCAmr.*` | future `src/amr/` | Not started |
| `Source/EB.*` | future `src/eb/` | Not started |
| `Source/LES.*` | future `src/les/` | Not started |
| `Source/Particle.cpp` | future `src/particles/` | Not started |

A row is called implemented only when its Fortran subsystem has an automated numerical gate. The current Riemann, characteristic, flattening, and transverse-correction rows remain qualified because general EOS, species, embedded boundaries, source terms, AMR, and 3D behavior are outside their tested scope.
