# PeleC to PeleF mapping

This table maps responsibilities, not source lines.

| PeleC reference | PeleF implementation | Status |
|---|---|---|
| `Source/main.cpp` | `app/pelef.F90`, `app/pelef2d.F90` | Separate serial 1D and 2D drivers |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base single-species state indices |
| PelePhysics constant-`gamma` calls | `src/physics/eos_ideal_mod.F90` | Existing hydro closure |
| PelePhysics species thermodynamics | `nasa7_thermo_mod`, `thermo_database_mod` | NASA7 H2/H/O/O2/OH/H2O/N2 subset verified |
| PelePhysics mixture EOS/caloric properties | `mixture_thermo_mod` | Independent ideal-gas mixture layer verified; hydro coupling pending |
| AMReX mesh/geometry responsibility | `mesh_mod`, `mesh_2d_mod` | Uniform 1D and Cartesian 2D meshes |
| AMReX boundary fill | `boundary_conditions_mod` and periodic wrapping in `ctu_2d_mod` | 1D outflow/periodic and 2D periodic subset |
| `Source/Riemann.H` LF path | `riemann_rusanov_mod` | Independent Rusanov implementation |
| `Source/Riemann.H` acoustic solver | `riemann_pelec_mod` | Qualified single-species ideal-gas subset |
| Direction-dependent flux assembly | `directional_flux_mod` | x/y momentum rotation and y fluxes verified |
| `Source/PLM.H` characteristic projection/tracing | `reconstruction_pelec_plm_mod` | Qualified 1D regular-cell subset |
| `Source/PLM.H::plm_slope` | `pelec_limited_slope` | Order 2 and 4 formulas verified |
| `Source/Godunov.H::flatten` | `pelec_flattening_coefficient` | 1D regular-cell formula verified |
| `Source/Godunov.*` transverse update responsibility | `src/hydro/ctu_2d_mod.F90` | Qualified periodic 2D regular-grid CTU-style subset |
| species conserved-state block | `multispecies_state_mod` | Passive runtime layout verified |
| passive/species Godunov fluxes | `multispecies_flux_mod`, `reconstruction_multispecies_mod` | Mass-flux closure and contact-wave tracing verified |
| multidimensional species update | `ctu_multispecies_2d_mod` | Periodic CTU subset verified |
| `Exec/RegTests/Sod` | `cases/sod`, `tools/compare_sod.py` | Exact-solution regressions |
| `Exec/RegTests/Shu-Osher` | `cases/shu_osher`, `tools/check_shu_osher.py` | Deterministic shock-wave gate |
| `Exec/RegTests/Sedov` | `cases/sedov`, `tools/check_sedov.py` | Independent planar strong-blast gate |
| `Exec/RegTests/MultiSpecSod` | `cases/multispec_sod`, `check_multispec_sod.py` | Passive multispecies regression implemented |
| `Source/React.cpp` reaction-source responsibility | `elementary_kinetics_mod`, `constant_volume_reactor_mod` | Elementary, third-body, Troe, and implicit 0D coupling verified |
| PelePhysics generated mechanism kernels | `generate_elementary_mechanism.py`, `src/generated/h2o2_*_mechanism_mod.F90` | Elementary and full normalized-JSON generation with cleanliness gates |
| reversible pressure-dependent chemistry | NASA7 equilibrium constants, third-body efficiencies, and Troe falloff | Full 29-reaction Cantera parity implemented |
| stiff reactor integration | dense backward Euler/Newton with adaptive step doubling | Verification implementation complete; CVODE/sparse path pending |
| third-body/falloff chemistry | `elementary_kinetics_mod` | Third-body and Troe falloff implemented; SRI pending |
| complete mechanism parsing | future Cantera YAML/CHEMKIN parser | Not started |
| `Source/PPM.*` | future `reconstruction_ppm_mod` | Not started |
| `Source/WENO.H` | future `reconstruction_weno_mod` | Not started |
| `Source/Diffusion.*` | future `src/diffusion/` | Not started |
| `Source/PeleCAmr.*` | future `src/amr/` | Not started |
| `Source/EB.*` | future `src/eb/` | Not started |
| `Source/LES.*` | future `src/les/` | Not started |
| `Source/Particle.cpp` | future `src/particles/` | Not started |

A row is called implemented only when its Fortran subsystem has an automated numerical gate. The current chemistry row is qualified: it covers the complete small Cantera H2/O2 mechanism and a dense constant-volume implicit solver, not a general PelePhysics parser, production-scale sparse integrator, hydrocarbon mechanism, or chemistry-coupled flow.
