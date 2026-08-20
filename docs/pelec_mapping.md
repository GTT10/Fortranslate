# PeleC to PeleF mapping

This table maps responsibilities, not source lines.

| PeleC reference | PeleF implementation | Status |
|---|---|---|
| `Source/main.cpp` | `app/pelef.F90`, `app/pelef2d.F90`, `app/pelef_reactive_1d.F90`, `app/pelef_reactive_2d.F90` | Separate constant-gamma and reactive serial drivers |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base single-species state indices |
| PelePhysics constant-`gamma` calls | `src/physics/eos_ideal_mod.F90` | Existing hydro closure |
| PelePhysics species thermodynamics | `nasa7_thermo_mod`, `thermo_database_mod` | NASA7 H2/H/O/O2/OH/H2O/N2 subset verified |
| PelePhysics mixture EOS/caloric properties | `mixture_thermo_mod`, `reactive_1d_mod`, `reactive_2d_mod` | NASA7 ideal-gas-mixture layer coupled to qualified reactive 1D/2D paths |
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
| reactive multidimensional Godunov/CTU responsibility | `reactive_2d_mod`, `pelef_reactive_2d` | Periodic regular-grid PCM/characteristic-PLM/characteristic-PPM normal predictors, full-state transverse correction, and x/y reduction verified |
| `Source/React.cpp` reaction-source responsibility | `elementary_kinetics_mod`, `constant_volume_reactor_mod`, `reactive_1d_mod`, `reactive_2d_mod` | Four-reaction 0D chemistry and Strang-split 1D/2D cell coupling verified |
| PelePhysics generated mechanism kernels | `generate_elementary_mechanism.py`, `src/generated/h2o2_elementary_mechanism_mod.F90` | Normalized JSON generation and cleanliness gate implemented |
| reversible elementary chemistry | NASA7 equilibrium constants and generated H2/O2 rates | Four-reaction Cantera parity implemented |
| reactive hydro state/flux path | `reactive_1d_mod`, `reactive_2d_mod`, `pelef_reactive_1d`, `pelef_reactive_2d` | NASA7 conversion, directional Rusanov/HLLC, PLM/PPM normal prediction, CTU, and Strang splitting verified |
| stiff reactor integration | future CVODE/SUNDIALS layer | Not started |
| third-body/falloff chemistry | future kinetics extensions | Not started |
| complete mechanism parsing | future Cantera YAML/CHEMKIN parser | Not started |
| `Source/PPM.*` regular-cell normal predictor | `reactive_1d_mod`, `reactive_2d_mod` characteristic PPM paths | Five-point reconstruction and `u-c/u/u+c` profile integration verified in 1D and as x/y normal predictors before 2D CTU correction |
| `Source/WENO.H` | future `reconstruction_weno_mod` | Not started |
| PelePhysics `Source/Transport/Simple.H` | `transport_database_mod`, `mixture_transport_mod`, `pelef_transport_probe` | Qualified dilute ideal-gas subset: Chapman--Enskog/Wilke/Mathur/mixture-averaged diffusion |
| PeleC `Source/Diffterm.H`, `Source/Diffusion.cpp` | `reactive_diffusive_flux_x`, `advance_reactive_transport` | Periodic 1D viscous, conductive, barodiffusive, correction-velocity, and enthalpy-flux subset verified |
| `Source/Diffusion.*` multidimensional/EB responsibility | future 2D/AMR diffusion modules | Not started |
| `Source/PeleCAmr.*` | future `src/amr/` | Not started |
| `Source/EB.*` | future `src/eb/` | Not started |
| `Source/LES.*` | future `src/les/` | Not started |
| `Source/Particle.cpp` | future `src/particles/` | Not started |

A row is called implemented only when its Fortran subsystem has an automated numerical gate. The current chemistry row is qualified: it covers reversible elementary reactions and one small constant-volume subset, not a complete PelePhysics mechanism, stiff integration, full transport parity, or multidimensional diffusive flow.

| PeleC responsibility | PeleF current mapping | Qualification |
|---|---|---|
| regular-cell PPM edge interpolation and monotonicity | `reactive_1d_mod` and `reactive_2d_mod` characteristic PPM paths | general-EOS mixture normal predictors in 1D and both 2D coordinate directions |
| `PPM.cpp` normal characteristic tracing | `reactive_1d_mod` and `reactive_2d_mod` `characteristic_ppm` paths | density/normal velocity/pressure frozen-composition projection; species and transverse velocities on the middle wave; 2D uses direction rotation |
| `Godunov.H::flatten` | `reactive_ppm_flattening_coefficient` | One-dimensional regular-cell pressure/compression detector verified |
| Colella--Woodward contact steepening | `reactive_ppm_contact_steepening_factor` | Separate bounded density/species subset; not a claim that current PeleC enables this option |


| PeleC/PelePhysics transport path | PeleF 0.17.0 |
|---|---|
| `Diffusion.cpp` coefficient/flux/divergence workflow | `reactive_transport_2d_mod` |
| `Diffterm.H` stress, Fourier, species enthalpy flux | directional face-flux kernels |
