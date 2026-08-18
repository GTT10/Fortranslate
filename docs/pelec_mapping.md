# PeleC to PeleF mapping

This table maps responsibilities, not source lines.

| PeleC reference | PeleF implementation | Status |
|---|---|---|
| `Source/main.cpp` | `app/pelef.F90`, `app/pelef2d.F90` | Separate serial 1D and 2D drivers |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base single-species state indices |
| PelePhysics constant-`gamma` calls | `src/physics/eos_ideal_mod.F90` | Existing hydro closure |
| PelePhysics species thermodynamics | `nasa7_thermo_mod`, `thermo_database_mod` | NASA7 H2/O2/H2O/N2 subset verified |
| PelePhysics mixture EOS/caloric properties | `mixture_thermo_mod` | Independent ideal-gas mixture layer verified; hydro coupling pending |
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
| `Exec/RegTests/MultiSpecSod` | `cases/multispec_sod`, `check_multispec_sod.py` | Passive multispecies regression implemented |
| `Source/PPM.*` | future `reconstruction_ppm_mod` | Not started |
| `Source/WENO.H` | future `reconstruction_weno_mod` | Not started |
| `Source/Diffusion.*` | future `src/diffusion/` | Not started |
| `Source/React.cpp` reactor-source responsibility | `isomerization_reactor_mod`, `app/pelef0d.F90` | Toy constant-volume scaffold only; detailed parity not claimed |
| mechanism parsing/code generation | future `src/chemistry/` and `tools/` | Not started |
| `Source/PeleCAmr.*` | future `src/amr/` | Not started |
| `Source/EB.*` | future `src/eb/` | Not started |
| `Source/LES.*` | future `src/les/` | Not started |
| `Source/Particle.cpp` | future `src/particles/` | Not started |

A row is called implemented only when its Fortran subsystem has an automated numerical gate. The current Riemann, characteristic, flattening, and transverse-correction rows remain qualified because general EOS, species, embedded boundaries, source terms, AMR, and 3D behavior are outside their tested scope.

## Multispecies responsibility mapping

| PeleC responsibility | PeleF implementation | Status |
|---|---|---|
| species conserved-state block | `src/core/multispecies_state_mod.F90` | passive runtime layout verified |
| passive/species Godunov fluxes | `src/hydro/multispecies_flux_mod.F90` | mass-flux closure verified |
| species PLM tracing | `src/hydro/reconstruction_multispecies_mod.F90` | 1D contact-wave subset verified |
| multidimensional species update | `src/hydro/ctu_multispecies_2d_mod.F90` | periodic CTU subset verified |
| `Exec/RegTests/MultiSpecSod` responsibility | `cases/multispec_sod`, `tools/check_multispec_sod.py` | independent regression implemented |
