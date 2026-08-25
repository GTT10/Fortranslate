# PeleC to PeleF mapping

This table maps responsibilities, not source lines.

| PeleC reference | PeleF implementation | Status |
|---|---|---|
| `Source/main.cpp` | serial `app/pelef*.F90` and distributed `app/pelef_mpi*.F90` drivers | Separate constant-gamma, reactive, and MPI verification applications |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base single-species state indices |
| PelePhysics constant-`gamma` calls | `src/physics/eos_ideal_mod.F90` | Existing hydro closure |
| PelePhysics species thermodynamics | `nasa7_thermo_mod`, `thermo_database_mod` | NASA7 H2/H/O/O2/OH/H2O/N2 subset verified |
| PelePhysics mixture EOS/caloric properties | `mixture_thermo_mod`, `reactive_1d_mod`, `reactive_2d_mod` | NASA7 ideal-gas-mixture layer coupled to qualified reactive 1D/2D paths |
| AMReX mesh/geometry responsibility | `mesh_mod`, `mesh_2d_mod` | Uniform 1D and Cartesian 2D meshes |
| AMReX boundary fill | `boundary_conditions_mod` and periodic wrapping in `ctu_2d_mod` | 1D outflow/periodic and 2D periodic subset |
| `Source/Riemann.H` LF path | `riemann_rusanov_mod` | Independent Rusanov implementation |
| `Source/Riemann.H` acoustic solver | `riemann_pelec_mod`, `reactive_1d_mod` | Constant-`gamma` subset plus qualified NASA7 mixture star-state/wave-interpolation path |
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
| `Source/React.cpp` reaction-source responsibility | `elementary_kinetics_mod`, `constant_volume_reactor_mod`, `reactive_1d_mod`, `reactive_2d_mod`, `mpi_reactive_1d_mod` | Elementary/full chemistry and serial/distributed Strang coupling verified |
| PelePhysics generated mechanism kernels | `generate_elementary_mechanism.py`, `src/generated/h2o2_elementary_mechanism_mod.F90` | Normalized JSON generation and cleanliness gate implemented |
| reversible elementary chemistry | NASA7 equilibrium constants and generated H2/O2 rates | Four-reaction Cantera parity implemented |
| reactive hydro state/flux path | `reactive_1d_mod`, `reactive_2d_mod`, `pelef_reactive_1d`, `pelef_reactive_2d` | NASA7 conversion, directional Rusanov/HLLC/PeleC-style fluxes, PLM/PPM normal prediction, CTU, and Strang splitting verified |
| stiff reactor integration | `constant_volume_reactor_mod` adaptive implicit backward-Euler path | Verified with generated Jacobian, step doubling, and rollback; not CVODE parity |
| third-body/falloff chemistry | `elementary_kinetics_mod`, `h2o2_full_mechanism_mod` | Third-body efficiencies, pressure falloff, and Troe verified |
| complete mechanism parsing | future Cantera YAML/CHEMKIN parser | Not started |
| `Source/PPM.*` regular-cell normal predictor | `reactive_1d_mod`, `reactive_2d_mod` characteristic PPM paths | Five-point reconstruction and `u-c/u/u+c` profile integration verified in 1D and as x/y normal predictors before 2D CTU correction |
| `Source/WENO.H` | `reconstruction_weno_mod`, multilevel `reactive_1d_mod` path | WENO5-JS, WENO5-Z, WENO7-Z, and WENO3-Z implemented as optional characteristic-PPM edge reconstruction; fixed formula-parity points and three-level AMR gates verified |
| PelePhysics `Source/Transport/Simple.H` | `transport_database_mod`, `mixture_transport_mod`, `pelef_transport_probe` | Qualified dilute ideal-gas subset: Chapman--Enskog/Wilke/Mathur/mixture-averaged diffusion |
| PeleC `Source/Diffterm.H`, `Source/Diffusion.cpp` | `reactive_diffusive_flux_x`, `advance_reactive_transport` | Periodic 1D viscous, conductive, barodiffusive, correction-velocity, and enthalpy-flux subset verified |
| `Source/Diffusion.*` multidimensional/AMR/EB responsibility | `reactive_transport_2d_mod` for regular 2D cells, `eb_reactive_transport_2d_mod` for EB levels, `amr_eb_transport_2d_mod`, `amr_eb_multilevel_transport_2d_mod`, and `amr_eb_multipatch_transport_2d_mod` for EB AMR synchronization, and `amr_reactive_1d_mod` for two-level 1D | Regular-grid 2D, adiabatic-slip EB, conservative single-patch, three-level, and separated sibling-patch EB AMR, and conservative two-level 1D AMR transport verified |
| AMReX distributed-box responsibility | `mpi_domain_1d_mod`, `mpi_reactive_transport_1d_mod`, `mpi_reactive_1d_mod`, `mpi_amr_patch_1d_mod`, `mpi_amr_sparse_patch_1d_mod`, `mpi_amr_eb_patch_2d_mod` | Uneven 1D blocks, subcycle-weighted 1D/EB AMR owner maps, rank-local 1D storage and sparse physics, direct 1D owner migration and traffic, plus 2D EB owner chemistry/hydro/transport/full-physics transactions, sparse numerical storage and average-down, sparse-input/output full-physics composition, and direct owner-local fine hydro verified for 1/2/4/8 ranks |
| `Source/PeleCAmr.*` hierarchy/synchronization subset | `amr_hierarchy_1d_mod`, `amr_multipatch_1d_mod`, `amr_patch_tree_1d_mod`, `amr_patch_tree_reactive_1d_mod`, `amr_regrid_1d_mod`, `amr_reactive_1d_mod`, `amr_multilevel_reactive_1d_mod`, `amr_multipatch_reactive_1d_mod`, `amr_eb_hierarchy_2d_mod`, `amr_eb_multilevel_2d_mod`, `amr_eb_multilevel_reactive_2d_mod`, `amr_eb_flux_register_2d_mod`, `amr_eb_reactive_2d_mod`, `amr_eb_regrid_2d_mod`, `reactive_eb_amr_2d_driver_mod`, `mpi_amr_patch_1d_mod`, `mpi_amr_sparse_patch_1d_mod` | Runnable solution-driven arbitrary-depth 1D branching, adjacent sibling exchange, complete sparse physics, owner-local tag planning, direct transactional topology migration, input-driven temperature-tagged two-level EB single- or multipatch lifecycles, and a public three-level EB hydro/chemistry time loop with static or dynamic-topology checkpointing and tag-driven finest movement across regular or EB-cut interfaces |
| `Source/EB.*`, `pc_umdrv_eb`, and AMReX-Hydro redistribution | `eb_geometry_2d_mod`, `eb_reactive_wall_flux_2d_mod`, `eb_reactive_redistribution_2d_mod`, `eb_reactive_reconstruction_2d_mod`, `eb_reactive_hydro_2d_mod`, `eb_reactive_transport_2d_mod`, `amr_eb_transport_2d_mod`, `amr_eb_multilevel_transport_2d_mod`, `amr_eb_multipatch_transport_2d_mod`, `reactive_eb_2d_driver_mod`, `amr_eb_hierarchy_2d_mod`, `amr_eb_flux_register_2d_mod`, `amr_eb_reactive_2d_mod`, `amr_eb_regrid_2d_mod`, `reactive_eb_amr_2d_driver_mod`, `pelef_reactive_eb_2d`, `pelef_reactive_eb_amr_2d` | Nodal geometry with cell/face centroids, reactive slip-wall pressure flux, PCM or active-stencil characteristic-PLM face-center Godunov fluxes, tangential interpolation to open-face centroids, FluxRedist, selectable zeroth- or second-order weighted StateRedist, adiabatic-slip molecular transport with single-patch, three-level, and separated sibling-patch diffusive reflux, active-cell chemistry splitting, and checkpoint-capable single-patch, separated multipatch two-level, or three-level EB AMR lifecycles with tag-driven finest movement; unsplit EB transverse prediction, fourth-order StateRedist, periodic/ghost neighborhoods, dynamic parent levels, non-outflow refined boundaries, and MPI remain |
| `Source/LES.*` | future `src/les/` | Not started |
| `Source/Particle.cpp` | future `src/particles/` | Not started |

A row is called implemented only when its Fortran subsystem has an automated
numerical gate. Chemistry includes a complete ten-species, 29-reaction H2/O2
path, but this is not a claim of arbitrary PelePhysics mechanism parsing,
hydrocarbon chemistry, CVODE parity, or full transport parity.

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


| PeleC boundary concept | PeleF 0.18.0 / 0.69.0 |
|---|---|
| physical ghost fill | `reactive_boundary_2d_mod` |
| impermeable wall pressure flux | `reactive_wall_flux_x/y` |
| wall Fourier/species flux | impermeable or prescribed zero-net-mass species flux with coupled enthalpy in `reactive_transport_2d_mod` |
| fixed inflow / extrapolated outflow | boundary primitive sampling |

| PeleC/PelePhysics chemistry concept | PeleF 0.19.0 |
|---|---|
| third-body and falloff rate evaluation | `elementary_kinetics_mod` |
| cell-local stiff reactor | `constant_volume_reactor_mod` implicit path |
| runtime mechanism selection | reactive 1D/2D application dispatch |

| Distributed responsibility | PeleF 0.24.0 |
|---|---|
| rank-local block ownership | `mpi_domain_1d_mod` uneven contiguous decomposition |
| periodic ghost fill | nonblocking state and temperature halo exchange |
| global timestep and diagnostics | communicator-wide min/max/sum reductions |
| ordered output | root `MPI_Gatherv` reconstruction |
| distributed reactive advance | `mpi_reactive_1d_mod` transactional Strang composition |
