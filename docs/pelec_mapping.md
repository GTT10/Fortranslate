# PeleC to PeleF mapping

This table maps responsibilities, not source lines.

All rows use PeleC commit
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` and the recursive revisions in
`references/pelec_baseline.json` unless a row or validation record names a
more specific reference.

| PeleC reference | PeleF implementation | Status |
|---|---|---|
| `Source/main.cpp` | serial `app/pelef*.F90` and distributed `app/pelef_mpi*.F90` drivers | Separate constant-gamma, reactive, and MPI verification applications |
| `Source/IndexDefines.H` | `src/core/state_indices_mod.F90` | Base single-species state indices |
| PelePhysics constant-`gamma` calls | `src/physics/eos_ideal_mod.F90` | Existing hydro closure |
| PelePhysics species thermodynamics | `nasa7_thermo_mod`, `thermo_database_mod` | NASA7 H2/H/O/O2/OH/H2O/N2 subset verified |
| PelePhysics mixture EOS/caloric properties | `mixture_thermo_mod`, `reactive_1d_mod`, `reactive_2d_mod`, `reactive_3d_mod` | NASA7 ideal-gas-mixture layer coupled to qualified reactive 1D/2D/3D paths |
| AMReX mesh/geometry responsibility | `mesh_mod`, `mesh_2d_mod`, `mesh_3d_mod` | Uniform 1D and Cartesian 2D/3D meshes |
| AMReX boundary fill | `boundary_conditions_mod`, periodic wrapping in `ctu_2d_mod`, and direct periodic face wrapping in `finite_volume_3d_mod` and `reactive_3d_mod` | 1D outflow/periodic and regular 2D/3D periodic subset |
| `Source/Riemann.H` LF path | `riemann_rusanov_mod` | Independent Rusanov implementation |
| `Source/Riemann.H` acoustic solver | `riemann_pelec_mod`, `reactive_1d_mod` | Constant-`gamma` subset plus qualified NASA7 mixture star-state/wave-interpolation path |
| Direction-dependent flux assembly | `directional_flux_mod` | x/y/z momentum rotations and directional fluxes verified |
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
| regular and static two-level 3D coordinate, directional-flux, Euler, chemistry, transport, and AMR responsibility | `mesh_3d_mod`, `directional_flux_mod`, `finite_volume_3d_mod`, `reactive_3d_mod`, `reactive_transport_3d_mod`, `amr_hierarchy_3d_mod`, `amr_reactive_3d_mod`, `amr_reactive_transport_3d_mod`, `mpi_amr_reactive_3d_mod`, `mpi_amr_sparse_reactive_3d_mod`, 3D configuration/problem modules, `pelef3d`, `pelef_reactive_3d`, `pelef_amr_reactive_3d`, `pelef_mpi_amr_reactive_3d` | Uniform periodic NASA7 multispecies PCM or characteristic-PLM SSPRK2 evolution with optional chemistry/transport, plus static single-patch AMR hydro, chemistry, and `r^2`-subcycled molecular transport with limited space/time coarse ghosts, six-face reflux/interface limiting, average-down, and transport-disabled schema-2/schema-3 plus transport-enabled schema-4/schema-5 restart; the public MPI path stores parent-aligned rank-local coarse/fine slabs and supports exact fixed/selected 1/2/4/8-rank parity and changed-rank transport continuation; dynamic topology, physical boundaries, scalable I/O, and 3D EB AMR remain open |
| `Source/EB.*` 3D geometry, wall, cut-cell divergence, stabilization, reaction, transport, and restart composition | `eb_geometry_3d_mod`, `eb_reactive_wall_flux_3d_mod`, `eb_reactive_redistribution_3d_mod`, `eb_reactive_hydro_3d_mod`, `eb_reactive_transport_3d_mod`, `reactive_eb_cfl_3d_mod`, `reactive_3d_mod`, `simulation_config_reactive_eb_3d_mod`, `reactive_eb_3d_driver_mod`, `reactive_eb_3d_checkpoint_mod`, `pelef_reactive_eb_3d` | Exact static x/y/z axis-plane metrics with cellwise discrete geometric identity, NASA7 slip-wall pressure balance, PCM zero-gradient hydro outer faces, common-residual transactional raw/FluxRedist/order-zero-StateRedist Euler publication, active-state hydro/transport bounds, optional masked elementary chemistry, StateRedist SSPRK2 molecular transport, transactional `R-T-H-T-R`, zero outer/embedded-wall transport flux, exact covered storage, H/O/N diagnostics, deterministic EB CSV, and versioned transactional serial restart with fixed schema 1 or selected bundle/policy/composition-bound schema 2; reaction-aware timestep control, general wall/FluxRedist transport, general geometry, high-order EB hydro, AMR, and MPI remain open |
| reactive multidimensional Godunov/CTU responsibility | `reactive_2d_mod`, `pelef_reactive_2d` | Periodic regular-grid PCM/characteristic-PLM/characteristic-PPM normal predictors, full-state transverse correction, and x/y reduction verified |
| `Source/React.cpp` reaction-source responsibility | `elementary_kinetics_mod`, `constant_volume_reactor_mod`, `reactive_1d_mod`, `reactive_2d_mod`, `mpi_reactive_1d_mod` | Elementary/full chemistry and serial/distributed Strang coupling verified |
| PelePhysics generated mechanism kernels | `ingest_cantera_mechanism.py`, `generate_elementary_mechanism.py`, `mechanisms/h2o2_cantera.yaml`, and generated H2/O2 modules | Pinned `ohmech` build-time Cantera YAML ingestion, normalized JSON/Fortran cleanliness, and generated thermo/kinetics/transport loading qualified; arbitrary runtime mechanisms remain open |
| PelePhysics generated mechanism compilation and selected execution | `cmake/MechanismBundle.cmake`, `CheckMechanismSource.cmake`, the selected 0D/regular-1D/2D/3D, AMR-1D, reactive-EB-2D/3D, static/dynamic single-patch two-level, dynamic two-level multipatch, static/dynamic three-level reactive-EB-AMR-2D, and regular-MPI-1D front-end templates, `selected_composition_mod`, `selected_mechanism_runtime_mod`, the shared 2D/3D/AMR/EB/MPI application drivers, `gas_transport_mod`, and `sundials_constant_volume_reactor_mod` | An optional normalized JSON bundle is regenerated into an isolated build-tree module, exercised through its declared ABI, advanced by an installed serial constant-volume application, and dispatched through installed regular-grid 1D, 2D, periodic 3D, serial reactive-1D AMR, serial single-level reactive-2D/3D EB, static/dynamic single-patch two-level, dynamic two-level multipatch, static/dynamic one-patch-per-level three-level reactive-2D EB AMR, and bounded regular-MPI-1D chemistry/transport applications with shared validation; native 0D and fixed/selected full-H2/O2 regular 1D/2D/3D, AMR-1D two-level/multilevel/multipatch, EB 2D plane-chemistry/circle-transport, EB 3D plane chemistry/transport, static/dynamic two-level, dynamic two-level multipatch, static/dynamic three-level EB-AMR outputs, and MPI 1/2/4-rank outputs retain byte-exact gates, selected planar EB 3D schema 2 and selected EB AMR 2D single-patch, multipatch, and three-level schemas bind restart to bundle/policy/composition, and optional pinned SUNDIALS 7.2.0 CVODE BDF adds rollback-safe reduced-state integration, Cantera trajectory/rate gates, and generation-safe sequential multi-context ownership; runtime parsing and selected MPI AMR/EB persistence remain open |
| reversible elementary chemistry | NASA7 equilibrium constants and generated H2/O2 rates | Four-reaction Cantera parity implemented |
| reactive hydro state/flux path | `reactive_1d_mod`, `reactive_2d_mod`, `reactive_3d_mod`, `pelef_reactive_1d`, `pelef_reactive_2d`, `pelef_reactive_3d` | NASA7 conversion, directional Rusanov/HLLC/PeleC-style fluxes, and elementary/full-H2O2 chemistry/transport splitting verified through periodic PCM or characteristic-PLM SSPRK2 3D; PPM and CTU remain qualified only in lower-dimensional paths |
| stiff reactor integration | `constant_volume_reactor_mod` adaptive implicit backward-Euler path and optional `sundials_constant_volume_reactor_mod` CVODE BDF path | Native step-doubling remains the default; optional official-Fortran SUNDIALS 7.2.0 uses a serial dense solver, largest-species reduced closure, energy-constrained semi-analytic Jacobian, caller rollback, deterministic fixture, pinned H2/O2 Cantera gates, and up to 64 generation-checked sequentially interleaved contexts; threaded/reentrant calls, sparse solvers, detailed fuels, and performance remain open |
| third-body/falloff chemistry | `elementary_kinetics_mod`, `h2o2_full_mechanism_mod` | Third-body efficiencies, pressure falloff, and Troe verified |
| complete mechanism parsing | bounded Cantera YAML importer plus future general Cantera/CHEMKIN path | One pinned ideal-gas NASA7/LJ elementary/three-body/Troe H2/O2 phase is qualified at build time; general YAML/CHEMKIN and runtime selection remain open |
| `Source/PPM.*` regular-cell normal predictor | `reactive_1d_mod`, `reactive_2d_mod` characteristic PPM paths | Five-point reconstruction and `u-c/u/u+c` profile integration verified in 1D and as x/y normal predictors before 2D CTU correction |
| `Source/WENO.H` | `reconstruction_weno_mod`, multilevel `reactive_1d_mod` path | WENO5-JS, WENO5-Z, WENO7-Z, and WENO3-Z implemented as optional characteristic-PPM edge reconstruction; fixed formula-parity points and three-level AMR gates verified |
| PelePhysics `Source/Transport/Simple.H` | `transport_database_mod`, `mixture_transport_mod`, `pelef_transport_probe` | Qualified dilute ideal-gas subset: Chapman--Enskog/Wilke/Mathur/mixture-averaged diffusion |
| PeleC `Source/Diffterm.H`, `Source/Diffusion.cpp` | `reactive_diffusive_flux_x`, `advance_reactive_transport` | Periodic 1D viscous, conductive, barodiffusive, correction-velocity, and enthalpy-flux subset verified |
| `Source/Diffusion.*` multidimensional/AMR/EB responsibility | `reactive_transport_2d_mod` and `reactive_transport_3d_mod` for regular cells, `amr_reactive_transport_3d_mod` and `mpi_amr_sparse_reactive_3d_mod` for static 3D AMR, `eb_reactive_transport_2d_mod` for 2D EB levels and isothermal/no-slip wall fluxes, `eb_reactive_transport_3d_mod` for planar 3D StateRedist transport, `amr_eb_transport_2d_mod`, `amr_eb_multilevel_transport_2d_mod`, and `amr_eb_multipatch_transport_2d_mod` for 2D EB AMR synchronization, and `amr_reactive_1d_mod` for two-level 1D | Periodic regular-grid 2D/3D transport, conservative serial/sparse-MPI static two-level 3D AMR transport, adiabatic-slip or first-order isothermal/no-slip 2D EB, adiabatic-slip/impermeable exact-axis-plane 3D EB StateRedist transport, conservative single-patch, three-level, separated sibling-patch and arbitrary-depth 2D EB AMR, sparse-MPI ownership, and conservative two-level 1D AMR transport verified; dynamic 3D AMR and 3D EB general-wall/AMR/distributed transport remain open |
| AMReX distributed-box responsibility | `mpi_domain_1d_mod`, `mpi_reactive_transport_1d_mod`, `mpi_reactive_1d_mod`, `mpi_amr_patch_1d_mod`, `mpi_amr_sparse_patch_1d_mod`, `mpi_amr_eb_patch_2d_mod`, `mpi_amr_reactive_3d_mod`, `mpi_amr_sparse_reactive_3d_mod` | Uneven 1D blocks, subcycle-weighted 1D/EB AMR owner maps, rank-local 1D storage and sparse physics, direct 1D owner migration and traffic, 2D EB owner full-physics transactions, and parent-aligned rank-local static 3D AMR state with two-layer hydro/transport halos, PCM/PLM hydro, `r^2` molecular transport, exact fixed/selected serial-rank parity, and rank-neutral root compatibility I/O; dynamic 3D topology and scalable I/O remain open |
| `Source/PeleCAmr.*` hierarchy/synchronization subset | `amr_hierarchy_1d_mod`, `amr_multipatch_1d_mod`, `amr_patch_tree_1d_mod`, `amr_patch_tree_reactive_1d_mod`, `amr_regrid_1d_mod`, `amr_reactive_1d_mod`, `amr_multilevel_reactive_1d_mod`, `amr_multipatch_reactive_1d_mod`, `amr_eb_hierarchy_2d_mod`, `amr_eb_patch_tree_2d_mod`, `amr_eb_patch_tree_reactive_2d_mod`, `amr_eb_multilevel_2d_mod`, `amr_eb_multilevel_reactive_2d_mod`, `amr_eb_flux_register_2d_mod`, `amr_eb_reactive_2d_mod`, `amr_eb_regrid_2d_mod`, `reactive_eb_amr_2d_driver_mod`, `mpi_amr_patch_1d_mod`, `mpi_amr_sparse_patch_1d_mod` | Runnable solution-driven arbitrary-depth 1D branching, adjacent sibling exchange, complete sparse physics, owner-local tag planning, direct transactional topology migration, input-driven temperature-tagged two-level EB single- or multipatch lifecycles, a public three-level EB hydro/chemistry time loop, and a separate arbitrary-depth branching 2D EB numerical tree with transactional migration, all-node CFL selection, recursive hydrodynamics and SSPRK2 transport, active-cell chemistry, atomic `R-H-R`/`R-T-H-T-R` splitting, and serial per-parent temperature-tag planning/rebuild |
| `Source/EB.*`, `pc_umdrv_eb`, and AMReX-Hydro redistribution | `eb_geometry_2d_mod`, `eb_reactive_wall_flux_2d_mod`, `eb_reactive_redistribution_2d_mod`, `eb_reactive_reconstruction_2d_mod`, `eb_reactive_hydro_2d_mod`, `eb_reactive_transport_2d_mod`, `amr_eb_transport_2d_mod`, `amr_eb_multilevel_transport_2d_mod`, `amr_eb_multipatch_transport_2d_mod`, `reactive_eb_2d_driver_mod`, `amr_eb_hierarchy_2d_mod`, `amr_eb_flux_register_2d_mod`, `amr_eb_reactive_2d_mod`, `amr_eb_regrid_2d_mod`, `reactive_eb_amr_2d_driver_mod`, `pelef_reactive_eb_2d`, `pelef_reactive_eb_amr_2d` | Nodal geometry with cell/face centroids, reactive slip-wall pressure flux, PCM or active-stencil characteristic-PLM face-center Godunov fluxes, tangential interpolation to open-face centroids, FluxRedist, selectable zeroth- or second-order weighted StateRedist, configurable adiabatic/isothermal slip/no-slip molecular transport with single-patch, three-level, separated sibling-patch, and arbitrary-depth diffusive reflux, active-cell chemistry splitting, restart-safe wall/transport fingerprints, selected-context-bound static/dynamic single-patch two-level, dynamic two-level multipatch, and static/dynamic one-patch-per-level three-level restart, and checkpoint-capable fixed-depth or sparse-MPI arbitrary-depth EB AMR lifecycles; unsplit EB transverse prediction, fourth-order StateRedist, periodic/ghost neighborhoods, multiple dynamic fixed-depth parents, non-outflow refined boundaries, and fixed-depth MPI distribution remain |
| `Source/LES.*` | future `src/les/` | Not started |
| `Source/Particle.cpp` | future `src/particles/` | Not started |

A row is called implemented only when its Fortran subsystem has an automated
numerical gate. Chemistry includes a complete ten-species, 29-reaction H2/O2
path, but this is not a claim of arbitrary PelePhysics mechanism parsing,
hydrocarbon chemistry, general CVODE mechanism coverage, or full transport
parity.

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

| Sparse MPI EB AMR responsibility | PeleF 0.200.0 |
|---|---|
| outflow-boundary recursive transport | four-level viscosity, Fourier conduction, mixture-averaged species diffusion, barodiffusion, correction velocity, and species enthalpy flux with physical-side register omission, serial/sparse rank parity, and cross-rank restart |
| outflow-boundary reacting full physics | owner-local elementary chemistry composed transactionally with recursive transport and hydro across fresh 1/2/4/8-rank runs and 2-to-4/8-rank restart |
| rank-local persistent state | root row tiles and exclusive fine-child payloads |
| coarse/fine restriction | targeted child-to-intersecting-root-owner buffers |
| root hydro and transport physics | owner-tiled finite-halo Euler stages retaining start, end, corrected state/temperature, and interface fluxes only on tile owners with no post-compute root assembly |
| fine-owner coarse context | direct patch-plus-two start/end/corrected state and temperature plus child-intersecting x/y flux fragments from root tile owners |
| coarse interface-flux consumption | globally indexed patch-local x/y face rectangles with a complete-root compatibility wrapper |
| coarse interface-flux ownership | retained root-tile x rows and unique y faces routed directly to the child owner for register accumulation |
| hydro child context and reflux | direct patch-plus-two start/end/corrected state and interface-flux fragments from tile owners, child-local context extraction and compact reflux, and direct corrected-fragment return |
| child exterior state-context extraction | globally indexed patch-plus-one start/end state and temperature support with a complete-root wrapper |
| reflux ordering | child-local support reflux, retained owner-local fine field, and ordered corrected fragments returned directly to intersecting root tile owners |
| final transport root commit | corrected tile state commits locally without root-owner row scatter |
| final hydro root commit | corrected tile state commits locally without root-owner row scatter |
| transport SSPRK2 blend | tile-local conserved-state average and EB-band EOS recovery with no root-field traffic |
| EB-cut conservation closure | tile-local physical-boundary flux contributions, communicator-wide conserved-vector sum, and tile-local correction |
| stable coarse timestep | owner-local EB hydro/transport limits, refinement scaling, and communicator minimum |
| public full-physics clock | repeated sparse stable-step selection, exact target clipping, and committed `R-T-H-T-R` accounting |
| explicit topology change | direct child-to-root restriction, distinct-new-owner PCM root assembly, overlap owner migration, and atomic one-copy commit |
| arbitrary-depth dynamic topology | per-parent owner-local temperature tags, domain-inclusive outflow-boundary children, compact plan reduction, caller EB geometry rebuild, deterministic redistribution, direct retained-overlap migration, and atomic commit |
| arbitrary-depth checkpoint/restart | selected-root gather/write, root-only read, compact geometry broadcast, domain-inclusive physical-side child reconstruction, owner-map recomputation, and direct rank-neutral scatter |
| arbitrary-depth composite output | one deterministic finest-available-cell CSV, with selected-root direct sparse gather and root-only file access |
| public serial arbitrary-depth application | dedicated namelist-driven root initialization/restart, recursive tag/regrid schedule, committed `R-T-H-T-R` clock, checkpoint calls, and composite output |
| public application restart parity | separate uninterrupted, checkpoint-stop, and restart processes with identity-keyed composite topology and field comparison |
| public branching application lifecycle | separated boundary-touching and interior tag features, multiple ordered patches on one level, four-level full physics, and serial/changed-rank restart parity |
| public sparse MPI arbitrary-depth application | namelist-driven sparse ownership, owner-local lifecycle, selected-root I/O, and 1/2/4/8-rank composite parity |
| public sparse MPI cross-rank restart | two-rank checkpoint-stop followed by independent four- and eight-rank restarts with ownership-weight changes and identity-keyed parity against an uninterrupted one-rank process |
| public sparse MPI fresh initialization | geometry-only topology and ownership first, reactive fields allocated on the sole root-node owner, then zero-copy allocatable transfer into sparse storage with no non-owner numerical root field |
| public checkpoint compatibility | schema-8 structured mesh, EB wall, transport physics, StateRedist, prolongation method, hierarchy, regrid fingerprint, cumulative transport-limiter diagnostic, original composite-integral baseline, fixed-capacity per-level chemistry/transport/hydro advance counters, and cumulative regrid-evaluation/tagged-cell history shared by serial and sparse-MPI applications with continuation and ownership controls explicitly excluded; fixed-depth formats store the same method in schema 3 |
| AMReX multilevel EB re-redistribution responsibility | topology-derived union of clipped three-by-three coarse/fine interface supports, excluding refined and covered parent cells, shared by fixed-depth, multipatch, arbitrary-depth, serial, and sparse-MPI conservation closures; exact AMReX transfer bookkeeping is not claimed |
| PeleC embedded diffusive wall flux | first-order centroid-normal isothermal Fourier flux and no-slip Newtonian traction/work in `eb_reactive_transport_2d_mod`, scaled by EB wall length and cut-cell fluid volume; exact quadratic-stencil parity is not claimed |
| public embedded-wall controls | shared single-level, fixed-depth AMR, arbitrary-depth AMR, and sparse-MPI `&embedded_boundary` kind, thermal mode, temperature, and velocity with explicit transport-dependency validation and restart fingerprinting |
| reactive EB AMR coarse-to-fine initialization | `pcm` or conservative limited `linear` selected by public input for fixed-depth and arbitrary-depth serial/sparse-MPI initialization, regrid, checkpoint, and restart; linear cut parents use connected 3-by-3 fluid-centroid least-squares slopes, volume-weighted zero-mean child offsets, component bounds, and EOS-admissibility PCM retry |
| checkpoint/output boundary | one packed payload per remote root tile or child gathered only to a selected root; non-root complete fields stay unallocated |
| formatted checkpoint and CSV output | selected root alone invokes the serial-compatible checkpoint writer and deterministic root/child CSV writers; completion status is collective |
| formatted checkpoint restart | selected root alone reads complete fields, then sends each root tile or child directly to its current sparse owner from a replicated geometry-only descriptor with no field broadcast or non-root child-field template |
| arbitrary-depth EB numerical tree | runtime relation sequence with ordered parent/child offsets, branching, per-level refinement ratios, transactional state migration, composite integration and CSV output, all-node hydro/transport stable-step reduction, recursive ratio-subcycled hydrodynamics and SSPRK2 molecular transport, per-child reflux, subtree conservation closure, active-cell chemistry, atomic `R-H-R`/`R-T-H-T-R` splitting, a stop-time-clipped committed-step clock, serial and owner-local MPI per-parent temperature-tag planning/rebuild, deterministic subcycle-weighted MPI node ownership, owner-only field allocation, direct owner migration and physics routing, and rank-neutral selected-root checkpoint/restart |

| Selected sparse MPI AMR 1D responsibility | PeleF 0.237.0 qualification |
|---|---|
| configure-time mechanism dispatch | selected NASA7, elementary kinetics, gas transport, exact-name bundle-order composition, and generated explicit/implicit policy enter one shared fixed/selected sparse MPI lifecycle |
| rank consensus | every rank must present the same bundle SHA-256, integrator, and normalized composition before initialization or file creation |
| chemistry policy | the generated integrator is forwarded through both sparse owner-local chemistry half-steps at every AMR level |
| selected restart provenance | schema 2 binds bundle SHA-256, integrator, species order, normalized composition, and original composite-integral baseline |
| changed-rank restart | checkpoint stores no communicator size or owner map; one-rank stop resumes on two and four ranks with exact uninterrupted output |

| Selected sparse MPI reactive EB patch-tree 2D responsibility | PeleF 0.238.0 qualification |
|---|---|
| shared application lifecycle | fixed and selected front ends enter `mpi_reactive_eb_patch_tree_2d_application_mod` for identical ownership, dynamic regrid, physics, diagnostics, and output control |
| selected initialization | exact-name bundle-order composition reaches the sole root owner and reactive boundary construction |
| rank consensus | bundle SHA-256, generated integrator, species/reaction/transport counts, and normalized composition agree before numerical field initialization |
| chemistry policy | the generated policy reaches both owner-local chemistry half-steps on every active patch-tree level |
| dynamic rank parity | fixed one rank and selected one, two, and four ranks produce byte-identical four-level composite CSV output |
| persistence boundary | fixed schema 8 is byte-frozen; every selected checkpoint/restart schedule or path is rejected pending a context-bound schema |

| Selected sparse MPI reactive EB patch-tree 2D persistence responsibility | PeleF 0.239.0 qualification |
|---|---|
| selected provenance | exclusive schema 9 binds bundle SHA-256, generated integrator, species order, normalized bundle-order composition, complete numerical fingerprint, and original composite baseline |
| transactional read | root validates context, geometry, topology, clock, counters, regrid history, physical fields, terminal marker, and strict end-of-stream before distribution |
| rank-neutral restart | checkpoint contains no communicator size or owner map; validated topology is redistributed and fields are scattered to owners for the active communicator |
| changed-rank parity | one-rank checkpoint-stop resumes on two and four ranks to exact uninterrupted fixed and selected final CSV bytes |
| I/O boundary | collective required-metadata presence is checked before root I/O or gather; root-formatted schema 9 is validation I/O, bundle SHA is provenance only, and crash-atomic replacement, payload authentication, and scalable I/O remain open |

| Rank-local static MPI AMR 3D responsibility | PeleF 0.240.0 qualification |
|---|---|
| coarse ownership | deterministic uneven contiguous x slabs with two periodic ghost planes |
| fine ownership | fine x planes aligned to parent coarse ownership; zero-fine-cell ranks are valid |
| persistent numerical storage | local coarse/fine conserved state and temperature only during stepping |
| high-order hydro | owner-local PCM or characteristic PLM SSPRK2 with MC/minmod limiting and time-interpolated limited coarse/fine ghosts |
| synchronization | deterministic six-face fine flux accumulation, local reflux, average-down, and EOS temperature recovery before collective commit |
| collective contract | ownership arrays, patch, NASA7 records, floating controls, solver, reconstruction, and limiter agree before payload-dependent communication |
| restart and output | unchanged root-formatted schema 2 with rank-neutral gather/scatter and exact two-to-four/eight-rank continuation; scalable I/O is not claimed |

| Fixed chemistry in static AMR 3D responsibility | PeleF 0.241.0 qualification |
|---|---|
| source composition | private coarse/fine cell reactors followed by average-down and covered-coarse temperature recovery |
| complete step | transactional `R(dt/2)-H(dt)-R(dt/2)` with exact disabled hydro dispatch |
| sparse ownership | owner-local coarse and fine source solves; ranks with no fine cells remain collective participants |
| collective context | complete ordered reactions, NASA7, patch/ownership, tolerances, integrator, solver, reconstruction, limiter, timestep, and chemistry-enable flag |
| invariants | five Euler composite integrals plus H/O/N totals; individual species must respond |
| parity | exact full-H2/O2 serial and one-/two-/four-/eight-rank coarse/fine CSV bytes |
| persistence boundary | chemistry-enabled schema-2 checkpoint/restart rejected; selected mechanisms and context-bound chemistry persistence remain open |

| Selected chemistry in static AMR 3D responsibility | PeleF 0.242.0 qualification |
|---|---|
| shared lifecycle | fixed and generated serial/MPI frontends enter mechanism-independent static two-level 3D AMR drivers |
| selected initialization | exact-name bundle-order composition, generated integrator, and selected NASA7/reaction/transport records reach the existing transactional `R-H-R` hierarchy step |
| pre-initialization MPI consensus | optional presence, bundle SHA-256, integrator, composition, ordered NASA7, reaction, and gas-transport records agree by exact character or real representation before state allocation or file creation |
| sparse participation | ranks with no fine planes remain in provenance, chemistry, CFL, advance, and acceptance collectives |
| diagnostic boundary | H/O/N totals are reported only for the qualified H2/O2/N2/AR species set; other selected mechanisms retain physical, closure, Euler, and synchronization checks |
| parity | the generated two-species fixture is exact across serial and one/two/four ranks; generated full-H2/O2 is exact against fixed serial and one/two/four/eight ranks |
| persistence and transport boundary | selected molecular transport, checkpoint, and restart requests are rejected before output; schema 2 remains fixed hydro-only compatibility I/O |

| Selected restart in static AMR 3D responsibility | PeleF 0.243.0 qualification |
|---|---|
| schema isolation | fixed calls retain schema 2; complete selected-context calls require schema 3 and cross-schema reads fail transactionally |
| provenance | bundle SHA-256, generated integrator, bundle-order composition, complete ordered reactions, chemistry tolerances, and NASA7 records |
| restored contract | solver/reconstruction/limiter, mesh/patch/CFL, physics flags, original composite baseline, clock/reflux history, and both conserved/temperature levels |
| transaction | private read candidates, EOS recovery, coarse/fine synchronization, terminal marker, and strict end-of-stream precede publication |
| MPI continuation | exact pre-I/O context consensus, root-formatted rank-neutral file, and two-rank checkpoint resumed exactly on one/two/four/eight ranks |
| I/O boundary | bundle digest is provenance only; transport, payload authentication, crash-atomic replacement, and scalable I/O remain open |

| Molecular transport in static AMR 3D responsibility | PeleF 0.244.0 qualification |
|---|---|
| transport physics | shared 3D Newtonian viscosity, Fourier conduction, mixture-averaged species diffusion, barodiffusion, correction velocity, and species-enthalpy flux |
| parabolic scheduling | hierarchy step `min(dt_c, r^2 dt_f)`; each SSPRK2 Euler stage advances the fine patch through `r^2` subcycles |
| fine boundary state | coarse start/end time interpolation with limited PLM coarse-to-fine transport ghosts on all six faces |
| synchronization | time- and area-averaged six-face fluxes, common uncovered-coarse species limiter, conservative fine-side flux correction, reflux, average-down, and EOS recovery |
| full physics | hierarchy-transactional `R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)` with exact rollback |
| sparse MPI | parent-aligned owner-local state, deterministic flux/theta reductions, complete transport-record and flag consensus, and valid zero-fine-plane ranks |
| parity | within each build type, fixed and selected full-H2/O2 output exact in serial and at one/two/four/eight ranks; active six-face limiter cases compare serial and MPI arithmetic directly |
| persistence boundary | transport-enabled checkpoint/restart rejected until a schema binds the gas-transport database and parabolic operator contract |

| Transport checkpoint persistence in static AMR 3D responsibility | PeleF 0.245.0 focused qualification |
|---|---|
| schema isolation | fixed public `full_h2o2` transport uses schema 4; selected chemistry transport uses schema 5; transport-disabled schema 2/3 remain compatible and cross-schema fallback is rejected |
| transport identity | complete ordered gas-transport database, parameter convention, operator identity, all controls, and cumulative transport diagnostics are serialized and validated |
| checkpoint boundary | writes occur at `POST_ACCEPTED_COARSE_STEP` after committed `R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)` stepping |
| transaction | private candidate read, terminal marker, strict EOF, and close precede publication; failed reads leave caller state unchanged |
| MPI continuation | all-rank context consensus precedes payload/root I/O; rank-neutral root files restart exactly on changed rank counts, including zero-fine-plane ranks |
| focused evidence | serial GNU Debug/Release `10/10` each and GNU/OpenMPI Debug/Release MPI `23/23` each; schema-4/5 SHA-256 values are in [`validation/0.245.0.md`](validation/0.245.0.md) |
| boundary | no 0.245 full matrix or clean-install result is claimed; crash-atomic replacement, payload authentication, scalable I/O, dynamic topology, physical boundaries, and 3D EB AMR remain open |
