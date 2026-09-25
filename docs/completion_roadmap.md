# Completion roadmap

## Definition and frozen reference

A selected PeleC responsibility is complete only when it has an independent
Fortran implementation and automated numerical evidence, or is explicitly
excluded with a reason and a usable alternative. A scaffold, test count or
clean compilation is not sufficient.

The reference is PeleC `development`
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; recursive dependencies are pinned in
[the manifest](../references/pelec_baseline.json). Changing that baseline must
identify affected numerical signatures and acceptance criteria.

The pinned reference rejects non-regular EB cells for Smagorinsky, Dynamic,
WALE and Vreman LES. Cut-cell LES is therefore a future model extension, not a
missing reference operator; the implemented regular-grid AMR LES path remains
the usable alternative. Independent serial/MPI SSPRK2 LES level stepping,
recursive intervals and flux synchronization now reach the MPI spray application
([coupling evidence](validation/les-mpi-spray-subcycle-20260925.md)). This does not
establish parity with the reference time integrator or physical validation.
See the [dated source audit](validation/les-reference-scope-20260925.md).

Watson latent heat and the PeleMP pressure-dependent boiling calculation now
have independent library APIs with distinct reference and boiling temperatures
([kernel evidence](validation/watson-properties-20260925.md)). The optional
`pelemp-cap` mode now connects distinct liquid inputs, the vapor-state cap and
restart identity, including ordinary EOS enthalpy-derived saturation latent heat
([integration evidence](validation/pelemp-cap-integration-20260925.md)).
The constant-latent `surface-cap` remains available and `reject` remains the
default. The ordinary EOS vapor-pressure path uses
enthalpy-derived latent heat, so replacing every latent heat with Watson would
not match the reference. Isolated property and fixed/selected regular
conservation/restart checks are complete. EB tree input and changed-rank
continuation now also have [focused evidence](validation/pelemp-eb-tree-input-20260925.md).
MPI AMR execution of this option, broader transfer-integrator parity and physical
qualification remain separate. Shared gas-species mapping now accumulates
component vapor and sources conservatively
([mapping evidence](validation/shared-vapor-mapping-20260925.md)). Independent
liquid molecular-weight inputs now feed liquid mole fractions, vapor pressure,
boiling and restart identity
([integration evidence](validation/liquid-molecular-weights-20260925.md)).
Single/mixed vapor-state conversion now uses the complete ambient mean molecular
weight and reference activity ordering
([skin-state evidence](validation/ambient-vapor-skin-20260925.md)).
GCM is excluded from this pinned PeleC port because upstream restricts it to
PeleLMeX; PeleMP is the available alternative. See the
[property source audit](validation/spray-property-scope-20260925.md).

## Consolidation versus solver completion

The source-provenance repair, documentation ownership, CMake module split,
archived PR11 supersession review, and hosted serial/MPI/CVODE/clean-install
workflows are implemented. [PR #102](https://github.com/GTT10/Fortranslate/pull/102)
records the actual integration status and completed run identities; see
[current status](current_status.md) and [the layout evidence](validation/build-layout-20260909.md).
A targeted NASA7 failure-status fix and its red/green tests are recorded
[separately](validation/nasa7-conversion-status-20260909.md).

Do not redo this housekeeping as a substitute for the remaining physical
capabilities. Conversely, do not label the following capabilities implemented
merely because the repository has been consolidated.

## 0.246.0 functional increment

The opt-in regular 3D path now connects coexisting single-component liquid populations selected per port, conservative
mass/momentum/total-energy exchange, drag/heating/evaporation, timed cone
injection from point or finite circular solid/hollow outlets with swirl with up to 16 independently configured ports, shared or independent constant/tabulated-rate schedules with prescribed or discharge-based speeds and optional per-port mass- or number-CDF Weibull/normal/lognormal/ChiSquared diameters, deviatoric Smagorinsky and coupled checkpoint/restart to the existing
gas solver. Cell-local CVODE is connected to this CFD path, with a pinned
FFCM1_Red liquid-methanol burn fixture and independent Cantera gas-only oracles.
Chemical labels, single-interval NASA7 and classified Lindemann/SRI/Tsang, chemically activated, PLOG and Chebyshev ingestion are
supported. Selected mechanism capacity is now 128 species; a 53-species GRI-Mech
0D and coupled-spray case has focused evidence. Large-mechanism MPI/scaling and
production diesel chemistry qualification remain open. Lu n-dodecane (53 species,
268 reactions) now imports with extended Troe coefficients and completes a
lean constant-volume ignition comparison; Newton trial-state handling preserves
accepted-state checks. See [focused evidence](validation/dodecane-20260915.md). See [scope and commands](spray.md) and [acceptance gates](validation/0.246.0.md).
These are bounded implemented capabilities. Multicomponent liquids and MPI/AMR
particle coupling are implemented; detailed-fuel physical qualification,
mixture boiling, general dispersion/collision and near-wall turbulence
qualification remain open.
Regular wall reflection and outflow now include gas/particle surface budgets
and restart persistence; isothermal/tangentially moving gas walls include molecular
and SGS heat/work budgets. The Ahamed impact-regime, splash-velocity, mass-distribution, four-child construction and deposited-film geometry kernels are connected by a transactional impact event with explicit wall impulse/energy budgets; regular coupled-step boundary handling now retains deposition patches, advances splash children and evolves stationary film temperature/mass with per-material constant or cubic conductivity; the CLI and coupled checkpoint retain these states, while physical qualification remains open. Lateral film transport is a future model extension beyond the pinned stationary-film path, not an omitted reference transport operator ([source audit](validation/film-reference-scope-20260918.md)). ETAB and KH-RT representative-drop
size/deformation, bounded computational parcel splitting and restart are implemented;
KH shed-mass accumulation and conservative four-child transverse release now
connect to CFD and restart, including computational splitting with partitioned
KH inventories and stationary wall impact processed before KH release on a
crossing substep. General turbulent dispersion, wet-impact correlation limits
and physical qualification remain open. Uniform gas/particle acceleration
now retains external impulse/work separately through coupled restart. Single-component
Antoine saturation-pressure fits are connected with explicit temperature bounds
and restart identity. Polynomial liquid density now drives injection, drag,
evaporation and breakup consistently; inverse-temperature liquid viscosity now
drives ETAB/KH-RT. The pinned piecewise drag law and ETAB deformation multiplier are selectable alongside Schiller-Naumann drag, with restart identity. An optional near-boiling vapor-state temperature cap follows
the pinned reference convention. Port material identity now routes evaporation to each liquid's gas species and persists through restart; density/pressure fits, vapor-state caps and breakup viscosity/surface tension can be set independently per port. Polynomial liquid heat capacity now supplies the thermal update and integrated liquid energy, including per-port fits. Mixed-liquid density, heat capacity, stored energy and fixed-composition heating
now have a shared thermodynamic API. Raoult surface partitioning and supplied
component-mass evaporation exchange, diffusion-weighted rates and a conservative
frozen-gas transport update are implemented. Regular-grid feedback, persistent
component masses, boundary budgets, representative breakup/splitting and composition
checkpoint APIs are connected. Combined hydro/transport/LES/chemistry stepping
now retains mixture composition through accepted updates and rollback. Physical mixture boiling rates remain open.
The mixture impact event kernel now partitions component masses into retained
parent, splash and film with wall budgets. Supplied film evaporation/heat exchange
now updates composition, temperature and height. Mixed film heat/evaporation rates
and a frozen-rate update are connected to grid feedback. Planar
crossing dispatch now computes mixed impact parameters and retains component budgets.
The checkpoint API retains mixed-film component masses. Stationary mixed-film grid
feedback now updates gas and film with wall budgets. Particle crossing, splash and
remaining-interval film evolution now share one transaction. Reactive stepping
retains both component matrices with full rollback. CLI wall inputs, component-film
output and restart are now connected. Shared-composition CLI inputs, initial/scheduled injection,
component output and restart are now connected. Mixed KH shed children preserve
composition and conservative transverse energy through coupling, splitting and restart.
Component density/heat-capacity fits, Antoine pressure and Clausius surface caps
are connected to CLI inputs and restart identity. Per-port compositions now feed
initial/scheduled injection, component budgets and density-based injection speeds.

## Ordered implementation work

EB hydro currently uses cell-centered Riemann states. The pinned PLM covered-cell
slope policy (4th/2nd/zero fallback) is now available through the normal PLM
reconstruction API ([checks](validation/eb-plm-slope-20260925.md)). An aperture-aware
three-direction normal predictor now returns primitive states and independently
traced internal energy ([checks](validation/eb-plm-normal-20260925.md)). A three-layer
halo API now constructs normal Riemann fluxes and Godunov pressure/velocity
([checks](validation/eb-plm-normal-flux-20260925.md)). Cell-local single/double
transverse correction now includes the joint aperture gate, pressure work and
primitive sources ([checks](validation/eb-transverse-20260925.md)). The twelve-solve
CTU flux pipeline is connected, including corrected Godunov pressure/velocity
([checks](validation/eb-ctu-flux-20260925.md)). A connectivity-gated face-centroid
flux interpolation API now follows the reference regular-area/aperture ordering
([checks](validation/eb-centroid-flux-20260925.md)). The prepared-halo hydro API
now assembles CTU, viscosity/normalization, centroid fluxes and conservative,
FluxRedist or StateRedist updates with integrated flux and inventory publication
([checks](validation/eb-ctu-update-20260925.md)). A cut-aware, multidirectional
flattening preparation API now feeds the prepared-halo update
([checks](validation/eb-flattening-20260925.md)). Physical wall/outflow Riemann
corrections and interior physical-wall viscosity masks now reach the complete
prepared CTU update ([checks](validation/eb-physical-flux-20260925.md)).
Configured plane/extruded/nodal geometry now supports domain-face extension,
and an outflow entry prepares state halos and flattening before the CTU update
([checks](validation/eb-extended-geometry-20260925.md)). Fixed/selected single-level applications now select `plm_ctu`, compose it with
reaction/transport and bind checkpoint continuation to the hydro method
([checks](validation/eb-ctu-application-20260925.md)). Serial/MPI patch-tree
halo APIs now accept six-layer target-resolution metrics with peer/ancestor
exchange and temporal interpolation ([checks](validation/eb-wide-halo-20260925.md)).
A conservative extended-geometry builder and serial common-interval CTU
patch stepping now connect peer synchronization and coarse/fine reflux
([checks](validation/eb-tree-ctu-20260925.md)). Serial recursive CTU stepping
and coupled reaction/transport composition now also connect time-interpolated
ancestor halos and accumulated fine fluxes
([checks](validation/eb-tree-ctu-subcycle-20260925.md)). MPI common CTU steps
now match serial under distributed and idle-rank ownership
([checks](validation/mpi-eb-ctu-common-20260925.md)). Recursive MPI CTU and
coupled common/recursive dispatch now also match serial
([checks](validation/mpi-eb-ctu-recursive-20260925.md)). Fixed/selected serial/MPI
tree application selection, scheduled regridding and changed-rank restart are
connected ([checks](validation/eb-tree-ctu-application-20260925.md)). Other physical
boundary halo preparation remains before claiming EB Godunov parity. MPI EB
spray now forwards CTU through coupled gas stepping and checkpoint identity
([checks](validation/mpi-eb-spray-ctu-20260925.md)). The serial single-level spray
entry and fixed/selected restart now also connect CTU
([checks](validation/serial-eb-spray-ctu-20260925.md)); broader physical-boundary
preparation and complete reference parity remain open. A single-level library
entry now prepares six-layer periodic/fixed-wall/prescribed-inflow halos for CTU
([checks](validation/eb-ctu-boundary-halos-20260925.md)); coupled transport,
application and tree/MPI boundary selection still need connection. The redistributed
transport Euler stage now accepts prepared donors and returns transactional
integrated face fluxes ([checks](validation/eb-transport-redistributed-halos-20260925.md));
SSPRK2 now invokes a caller-supplied donor provider at both stages and returns
weighted integrated fluxes ([checks](validation/eb-transport-stage-halos-20260925.md)).
Single-level coupled CTU now forwards inflow/outflow labels and a transport donor
provider with combined hydro/transport boundary accounting
([checks](validation/eb-coupled-inflow-boundaries-20260925.md)). Built-in inflow/outflow donors now reuse the regular boundary filler at every
transport stage ([checks](validation/eb-built-in-inflow-donors-20260925.md)).
Fixed/selected single-level gas applications now connect inflow/outflow input,
open-face donor timestep bounds and boundary-bound restart
([checks](validation/eb-inflow-application-20260925.md)). Fixed-wall transport policies
and isothermal heat flux now connect to single-level coupled CTU
([checks](validation/eb-fixed-wall-transport-20260925.md)). An aperture/fluid-volume
weighted isothermal timestep bound is available
([checks](validation/eb-isothermal-timestep-20260925.md)). Reflected configured
geometry and fixed/selected single-level wall input, thermal timestep and restart
binding are connected ([checks](validation/eb-wall-application-20260925.md)).
The periodic transport operator now wraps gradients and shares aperture-weighted
seam fluxes and species limiter coefficients
([checks](validation/eb-periodic-transport-20260925.md)). Periodic configured geometry,
coupled gas dispatch and fixed/selected application restart are now connected
([checks](validation/eb-periodic-application-20260925.md)). Fixed/selected single-level gas applications now accept independent inlet mixtures;
single-level spray state generation and restart now forward those mixtures. Serial single-level spray boundary
input now covers inflow/outflow, periodic and domain-wall gas conditions with
specular particle reflection ([checks](validation/eb-spray-wall-20260926.md)).
The serial EB tree halo API and common/recursive time updates now accept
prescribed inflow states. MPI gas stepping and fixed/selected application
dispatch now include changed-rank restart
([checks](validation/eb-tree-inflow-application-20260926.md)). Tree applications
also accept independent named mixtures per inlet face
([checks](validation/eb-tree-inlet-composition-20260926.md)). Periodic gas and MPI
spray boundaries are connected; adiabatic tree gas wall inputs and restart are connected
([checks](validation/eb-tree-wall-application-20260926.md)). Isothermal tree gas and MPI spray input, stepping and restart are connected
([checks](validation/mpi-eb-isothermal-spray-20260926.md)); MPI spray specular wall routing and changed-rank application
restart are connected
([checks](validation/mpi-eb-wall-spray-application-20260926.md)).

The 3D reactive flow path still lacks complete characteristic PPM CTU tracing.
A periodic three-axis normal PPM face-flux predictor now reuses the existing
1D/2D characteristic reconstruction, with analytic-flux, rotational and periodic
conservation [checks](validation/ppm-normal-3d-20260925.md). Twelve-solve
transverse/corner flux construction is now implemented with PeleF conservative
EOS recovery ([focused checks](validation/ppm-ctu-3d-20260925.md)); reference
primitive-pressure parity remains open. The periodic conservative CTU time
update and full-split library dispatch now have multistep conservation and
rollback [checks](validation/ppm-step-3d-20260925.md). Fixed/selected regular
gas application input now selects periodic PPM
([input checks](validation/ppm-input-3d-20260925.md)). Four-layer CTU halos and
physical-boundary library stepping now have exact subdomain-flux and boundary
ledger [checks](validation/ppm-halo-3d-20260925.md). Fixed/selected regular gas
applications now also select physical-boundary PPM with multistep molecular
transport [checks](validation/ppm-boundary-input-20260925.md). Serial common-step
AMR tree hydro now connects PPM and reflux
([checks](validation/ppm-amr-common-20260925.md)). Serial recursive level
stepping also supports PPM with integrated fluxes and 1/2/6-step synchronization
([checks](validation/ppm-amr-recursive-20260925.md)). MPI common-step PPM now
matches serial fields and diagnostics, including empty ranks
([checks](validation/ppm-mpi-amr-common-20260925.md)). Recursive MPI PPM now
also matches serial fields, diagnostics and 1/2/6 step counts
([checks](validation/ppm-mpi-amr-recursive-20260925.md)). Fixed/selected
serial/MPI AMR tree gas application input now enables PPM, with reacting
transport and fixed changed-rank continuation
[checks](validation/ppm-amr-application-20260925.md). MPI AMR spray input now also enables PPM,
with fixed 2-to-4-rank particle/gas continuation
[checks](validation/ppm-amr-spray-20260925.md). Regular serial spray also exposes PPM with Mach-boundary/transport/Smagorinsky
restart [checks](validation/ppm-regular-spray-20260925.md). Serial two-level AMR
hydro and reaction splitting now support PPM, including periodic/physical
boundaries, reflux and exact complete-refinement parity
([checks](validation/ppm-two-level-hydro-20260925.md)). Serial fixed/selected
two-level application input, conservative regrid, transport and checkpoint
selection now support PPM, with fixed reacting restart and exact selected output
[checks](validation/ppm-two-level-application-20260925.md). Sparse two-level MPI coarse PPM and four-layer halos now match serial even with
one owned cell per rank ([checks](validation/ppm-mpi-two-level-coarse-20260925.md)).
Fine-level MPI PPM now supports subcycling, reflux, reaction/transport splitting
and exact serial hierarchy parity
([checks](validation/ppm-mpi-two-level-hierarchy-20260925.md)). Fixed/selected MPI
two-level input and regridding also enable PPM; fixed four-to-two-rank restart
matches serial output ([checks](validation/ppm-mpi-two-level-application-20260925.md)).
A primitive-state Riemann entry now preserves traced pressure and uses the pinned
PeleC contact/regularization policy
([checks](validation/ppm-primitive-riemann-20260925.md)). An optional CTU API
now applies separate single/double transverse pressure work using Godunov face
pressure and velocity ([checks](validation/ppm-transverse-pressure-20260925.md)).
The pressure-traced path now uses the minimum shock flattening across all three
axes ([checks](validation/ppm-multidirectional-flattening-20260925.md)). Pinned
normal tracing recomputes internal energy from pressure through the EOS, as does
PeleF; independent normal internal-energy tracing is not an additional missing
reference step. The existing `characteristic_ppm`/`pelec` application combination
now selects this pressure policy through the common periodic/halo CTU entries
([checks](validation/ppm-pressure-dispatch-20260925.md)). Scalar reconstruction
and profile integration now agree with compiled pinned C++ bodies in a bounded
675-case comparison ([checks](validation/ppm-pinned-kernels-20260925.md)).
The primitive Riemann entry also agrees with the pinned compiled kernel using a
NASA7 EOS adapter in 56 cases ([checks](validation/riemann-pinned-kernel-20260925.md)).
Double-transverse pressure subtraction now follows the pinned spatial order
([checks](validation/ppm-transverse-order-20260925.md)). Regular-grid inflow
faces now receive final fourth-order-slope PLM flux replacement
([checks](validation/ppm-inflow-plm-20260925.md)). Serial/MPI two-level AMR now
forwards physical inflow masks on coarse/fine ownership boundaries and normalizes
PLM face composition before EOS evaluation
([checks](validation/ppm-two-level-inflow-20260925.md)). Serial/MPI patch-tree
common and recursive stepping now also forwards physical inflow masks
([checks](validation/ppm-tree-inflow-20260925.md)). Regular MPI spray hydro now
forwards domain-end inflow masks too
([checks](validation/ppm-regular-mpi-inflow-20260925.md)). External
normal primitive projection now preserves independent species weights and agrees
in a bounded 315-case compiled-kernel comparison
([checks](validation/ppm-mixture-projection-20260926.md)). Regular zero-source
single/double transverse kernels also agree in 144 compiled-reference cases
([checks](validation/ppm-pinned-transverse-20260926.md)). Godunov contact-speed
species upwinding and zero-viscosity final normalization now connect through
regular/EB CTU ([checks](validation/ppm-passive-flux-20260926.md)). The local
double-transverse kernel now accepts primitive sources and agrees in 216 compiled
reference cases ([checks](validation/ppm-transverse-source-20260926.md)); full
periodic/halo CTU source routing now also agrees in 28 full-assembly cases
([checks](validation/ppm-ctu-source-20260926.md)). Cell-local conserved-source conversion
is implemented ([checks](validation/ppm-source-conversion-20260926.md)). Grid
conversion now feeds periodic/halo CTU directly
([checks](validation/ppm-conserved-source-20260926.md)). Periodic CTU stepping
now adds supplied frozen conserved rates with inventory and rollback checks
([evidence](validation/ppm-source-step-20260926.md)). Boundary stepping also
accepts caller-supplied source halos and separates boundary/source inventories
([checks](validation/ppm-boundary-source-step-20260926.md)). The full library
split now forwards these sources with exact transport-composition checks
([evidence](validation/ppm-source-split-20260926.md)). A regular-cell external
force/work constructor now follows the pinned force-per-volume policy
([checks](validation/external-force-source-20260926.md)). Old/new external work
correction and a separate volume-source ledger now connect through full splitting
([checks](validation/external-force-step-20260926.md)). Fixed/selected regular
gas input and source-accounted diagnostics now connect external forcing
([checks](validation/external-force-app-20260926.md)); broader application routing
and full SDC integration remain open. The reference-ordered SDC candidate
assembler is now implemented and checked for source/hydro/reaction component
mapping ([checks](validation/sdc-state-assembly-20260926.md)). Nonreacting periodic
iteration now connects source callbacks, CTU and repeated state assembly
([checks](validation/sdc-periodic-iteration-20260926.md)). Fixed/selected regular
gas input selects this nonreacting periodic path with external force and source
ledgers ([checks](validation/sdc-app-20260926.md)). Molecular transport now
participates in old/new source evaluation with application inventory checks
([evidence](validation/sdc-transport-20260926.md)); forced chemistry and broader
boundary/AMR/MPI support remain open. The forced chemistry RHS now supports
variable density and externally changing internal energy
([checks](validation/forced-reactor-rhs-20260926.md)). An adaptive implicit
substep now integrates the forced reactor
([checks](validation/forced-reactor-integrator-20260926.md)). Full-interval
subcycling now preserves caller state on failure
([checks](validation/forced-reactor-interval-20260926.md)). Periodic SDC library
reaction coupling is implemented with a uniform reacting check
([evidence](validation/sdc-reactions-20260926.md)). Fixed/selected application
routing now retains reaction sources across steps
([checks](validation/sdc-reacting-app-20260926.md)); nonuniform reacting
qualification and broader boundary/AMR/MPI routing remain open. The library
now handles zero-gradient outflow with a separate boundary inventory
([checks](validation/sdc-outflow-20260926.md)). Fixed/selected dispatch now
supports reacting outflow with separate boundary/source accounting
([checks](validation/sdc-outflow-app-20260926.md)). Outflow molecular transport
now participates in old/new source evaluation and boundary accounting
([checks](validation/sdc-outflow-transport-20260926.md)); other physical
boundary types remain open.
SDC prescribed outlet pressure now reaches hydro and transport
([checks](validation/sdc-pressure-outlet-20260926.md)). Specified-state and Mach
inflow now share existing ghost and transport treatment in SDC
([checks](validation/sdc-inflow-20260926.md)). Slip, no-slip and isothermal
wall dispatch now reaches SDC hydro and transport
([checks](validation/sdc-walls-20260926.md)); AMR/MPI/EB SDC remain open. Physical
wall/outflow Riemann-state correction now reaches regular serial/MPI AMR and
MPI spray dispatch ([checks](validation/ppm-physical-riemann-20260926.md)).
The full periodic zero-source twelve-solve assembly now agrees in six bounded
compiled-reference cases ([checks](validation/ppm-pinned-ctu-20260926.md)).
The expanded ten-case oracle also covers active flattening and artificial
viscosity ([checks](validation/ppm-pinned-shock-20260926.md)). Eighteen cases now
include inflow/outflow/slip/no-slip boundary assembly; inflow PLM replacement
preserves independent species weights
([checks](validation/ppm-pinned-boundary-20260926.md)). Broader qualification
remains open.
The reference also adjusts final fluxes with compressive artificial viscosity
(default `difmag=0.1`). A regular-cell nodal-divergence/flux-correction API now
implements this step with no-slip tangential exclusions and species-flux
normalization ([checks](validation/artificial-viscosity-20260925.md)). Periodic and
halo CTU now apply coefficient 0.1 by default with the `pelec` pressure policy;
regular, AMR, tree and MPI spray callers provide physical no-slip masks
([checks](validation/ppm-viscosity-dispatch-20260925.md)). The CTU API accepts an
explicit coefficient. Fixed/selected regular gas namelist input now exposes it
([checks](validation/difmag-regular-input-20260925.md)); those executables have
no restart interface. Serial fixed/selected regular spray now forwards the
coefficient through coupled stepping and binds its resolved value to restart
([checks](validation/difmag-spray-input-20260925.md)). Regular MPI spray now forwards
collectively agreed coefficients and binds restart context
([checks](validation/difmag-mpi-spray-20260925.md)). Serial/MPI patch-tree common,
level and recursive hydro now forward explicit coefficients with collective
agreement and rollback ([checks](validation/difmag-amr-tree-hydro-20260925.md)).
Serial/MPI reaction/transport composition and the AMR spray gas call now forward
the coefficient ([checks](validation/difmag-amr-split-20260925.md)). Serial/MPI
tree gas and MPI AMR spray now expose coefficient input and bind
resolved restart identity ([checks](validation/difmag-amr-input-20260925.md)).
The older two-level AMR hydro, reaction and transport composition now also
forward explicit coefficients through coarse/fine stepping and reflux
([checks](validation/difmag-two-level-20260925.md)). Fixed/selected serial/MPI
two-level applications now expose input and bind the resolved coefficient to
restart identity ([checks](validation/difmag-two-level-input-20260925.md)).
Regular serial PCM/PLM library stepping now applies explicit positive coefficients
at both SSPRK2 stages, including mean face fluxes and boundary accounting
([checks](validation/pcm-plm-viscosity-20260925.md)). Explicit zero retains the unadjusted PCM/PLM calculation. Serial/MPI tree common and level hydro now apply
the coefficient before flux integration and reflux, including recursive
subcycling ([checks](validation/pcm-plm-tree-viscosity-20260925.md)). PCM/PLM
application defaults now resolve to the reference 0.1 through shared regular
gas configuration, including AMR checkpoint identity
([checks](validation/regular-gas-default-difmag-20260925.md)). Serial/MPI regular spray
now also defaults to 0.1 with changed-rank restart evidence
([checks](validation/spray-default-difmag-20260925.md)). EB application inputs also now default to 0.1 with zero-to-default restart
rejection ([checks](validation/eb-default-difmag-20260925.md)). Regular boundary/full-composition and PPM
flux APIs now default to 0.1 independently of solver/reconstruction
([checks](validation/regular-library-default-difmag-20260925.md)). Periodic
PCM/PLM direct stepping and the older Strang wrapper now accept the coefficient
and default to 0.1, preserving explicit-zero dispatch
([checks](validation/periodic-api-difmag-20260925.md)). Serial/MPI AMR tree
common/level hydro now resolves omission to 0.1, with common/recursive
serial-versus-MPI equivalence checks
([evidence](validation/amr-tree-default-difmag-20260925.md)). Serial/MPI
older two-level hydro now also resolves and forwards 0.1 on every coarse/fine
stage, with eight-rank default-versus-explicit checks
([evidence](validation/two-level-default-difmag-20260925.md)). Regular MPI
spray hydro/coupling library defaults now also resolve to 0.1, with four-rank
PLM/PPM omitted-versus-explicit checks
([evidence](validation/mpi-spray-library-default-20260925.md)). EB hydro and
hierarchy/coupling libraries now also default to 0.1; transport-only shared
operators retain zero
([evidence](validation/eb-library-default-20260925.md)). The viscosity kernel now accepts
EB volume fractions and open-face masks, suppressing nodal divergence near
cut/covered cells using the reference cutoff
([checks](validation/eb-viscosity-kernel-20260925.md)). The EB face-flux API and
plain/flux-redistributed/state-redistributed Euler APIs now accept explicit
coefficients. Outflow halos are built locally; external positive-coefficient
halos must include corner states and matching volume fractions
([checks](validation/eb-viscosity-dispatch-20260925.md)). Serial common-step EB
tree hydro now requests complete donor halos and forwards explicit coefficients
([checks](validation/eb-tree-viscosity-20260925.md)). Extracted patches now retain same-resolution exterior volume fractions through
serial/MPI halo fill and checkpoint restart
([checks](validation/eb-fine-geometry-halo-20260925.md)). Geometry without that
context still uses donor fractions; external reference geometry fidelity remains open,
as does the uniform default policy. MPI common-step hydro now
fetches complete owner-local halos and applies agreed coefficients, with
three-rank serial parity
([checks](validation/mpi-eb-tree-viscosity-20260925.md)). Serial recursive hydro and full
reaction/transport composition now forward the coefficient and request
time-interpolated complete halos
([checks](validation/eb-tree-viscosity-split-20260925.md)). MPI recursive and
coupled hydro now do likewise with collective coefficient agreement and
three-rank reacting/transport parity
([checks](validation/mpi-eb-viscosity-split-20260925.md)). Single-domain EB gas
input now accepts finite nonnegative coefficients and binds positive values
to checkpoint identity; omitted/zero retains legacy identity
([checks](validation/eb-difmag-input-20260925.md)). EB tree gas input now forwards
the coefficient through serial/MPI composition and patch checkpoint identity,
with a focused serial recursive restart check
([checks](validation/eb-tree-difmag-input-20260925.md)). Single-domain EB spray
now forwards explicit coefficients through particle/gas composition and binds
them to restart context while retaining the legacy zero component sequence
([checks](validation/eb-spray-difmag-20260925.md)). MPI EB spray now forwards
agreed coefficients before gas stepping, with 2/4-rank application parity
and changed-rank restart checks
([checks](validation/mpi-eb-spray-difmag-20260925.md)). Older two-level EB
common/recursive hydro now accepts explicit coefficients with complete patch
halos and conservative reflux
([checks](validation/eb-two-level-viscosity-20260925.md)). Its reaction/transport
split wrappers now forward the coefficient through common/recursive hydro;
see [coupled checks](validation/eb-two-level-viscosity-full-20260925.md).
Multipatch hydro/full composition now forwards coefficients through complete
peer halos, flux reconciliation and reflux
([checks](validation/eb-multipatch-viscosity-20260925.md)). The older
two-level/multipatch application now accepts explicit coefficients and binds
patch restart identity
([checks](validation/eb-two-level-difmag-input-20260925.md)). Regular MPI spray PLM hydro now
applies explicit coefficients with short-slab serial parity
([checks](validation/plm-mpi-spray-viscosity-20260925.md)); serial spray and MPI
mixture wrappers accept those coefficients. Serial two-level AMR now applies
explicit PCM/PLM coefficients on coarse and fine stages before flux averaging
and reflux ([checks](validation/pcm-plm-two-level-viscosity-20260925.md)). MPI
two-level AMR now applies the same corrections using local coarse/fine halos,
with eight-rank serial parity
([checks](validation/mpi-two-level-viscosity-20260925.md)). Shared gas configuration
now accepts explicit PCM/PLM coefficients, and two-level/tree checkpoints bind
them while preserving legacy zero identity
([checks](validation/nonppm-difmag-input-20260925.md)). Regular serial/MPI spray
now accepts explicit PLM coefficients and binds restart identity, with fixed
coupled and changed-rank checks
([evidence](validation/spray-plm-difmag-input-20260925.md)). Regular MPI spray
now exposes PPM through particle coupling and checkpoint context, with fixed
changed-rank continuation [checks](validation/ppm-regular-mpi-application-20260925.md). The regular MPI spray hydro library now accepts PPM
with four-layer exchange across two-cell slabs and serial state parity
[checks](validation/ppm-regular-mpi-hydro-20260925.md); the existing
PCM/PLM solver remains available.

| Order | Work unit | Measurable exit gate |
| --- | --- | --- |
| 1 | Validate existing physics before expanding it | Pin external cases, states, observables and tolerances; compare H2/O2 trajectories/rates with Cantera and flow fields with frozen PeleC; characterize convergence and conservation |
| 2 | Production chemistry path | Specify the target fuel mechanism and unsupported forms; validate ingestion/0D ignition; qualify stiff CFD integration with rollback and conservation before detailed-fuel claims |
| 3 | General 3D flow and mesh | Physical boundaries, then dynamic/multipatch/multilevel AMR and general EB in separate increments; require conservative regrid/reflux, changed-rank restart and convergence |
| 4 | LES and Lagrangian spray | Validate closures and particle submodels, then two-way mass/momentum/energy coupling; qualify evaporation and breakup before penetration/ignition comparisons |
| 5 | Practical release qualification | Reproducible clean install, documented supported inputs, distributed I/O/scaling measurements, physical benchmarks and owner-resolved licensing |

These are dependencies, not permission to bundle several physical models into
one PR. Existing regular 3D and static two-level AMR foundations must not be
reimplemented. Their restrictions remain explicit in current status.

Regular 3D hydro/chemistry/molecular transport now accepts per-face prescribed reservoir inflow, outflow, slip-wall and adiabatic/isothermal no-slip conditions with prescribed tangential wall velocity
with integrated boundary-flux accounting; see the [case guide](../cases/reactive_boundaries_3d/README.md).
Mach-state inflow is connected for regular and patch-tree gas/spray paths;
broader characteristic-boundary qualification and moving geometry remain open.
EB tree halos and serial common-step hydro now also accept Mach inlets
([focused evidence](validation/eb-mach-hydro-20260925.md)). Serial common/recursive
transport, recursive hydro and full splitting now forward those controls
([split evidence](validation/eb-mach-split-20260925.md)). MPI halos, common and
recursive hydro/transport and full splitting now also forward these controls
([MPI evidence](validation/eb-mach-mpi-20260925.md)). Serial/MPI timestep selection
now includes inlet signal speeds and diffusivity
([bound checks](validation/eb-mach-timestep-20260925.md)). Fixed/selected gas
application input and checkpoint binding are connected; fixed serial/2/3-rank
execution and rank-changing restart have
[focused evidence](validation/eb-mach-application-20260925.md). MPI EB spray now
forwards the same controls, with 2/4-rank conservation/parity and 2-to-4-rank
restart [checks](validation/eb-mach-spray-20260925.md). These bounded tests do
not establish nozzle validation or deposition/splash/film inlet qualification.
AMR-particle boundary coupling and stationary level-set EB applications are
implemented within the documented subsets. A three-axis extrusion API now maps
validated 2D level-set sections into 3D volume/face/wall metrics, including faceted
normal integrals and uniform-pressure balance. Nonplanar section-file geometry now connects
to fixed/selected application inputs and geometry-bound restart. Outward-neighbor StateRedist connects
nonplanar hydro/transport/reaction updates in a focused extruded-section case.
A 3D nodal level-set builder now supplies tetrahedrally reconstructed cell/face
metrics and pressure balance. Nodal-file geometry now connects to fixed/selected
application/restart and coupled reaction/transport/hydro updates. MPI EB spray
regridding and changed-rank restart are connected
([dated evidence](validation/spray-mpi-eb-regrid-application-20260918.md)); broader
small-cell stabilization, convergence and physical qualification remain open.

Serial/MPI 3D temperature-tagged patch movement/resizing is implemented with conservative
remapping, overlap reuse and geometry-aware changed-rank restart, including
physical-boundary tagging and boundary-touching patches. Remapping now runs
within fixed coarse slabs without full-level root aggregation. General
patch-tree evolution and integration are implemented as described below;
measured scaling and broader distributed I/O qualification remain open.
Two-level MPI AMR CSV output and final diagnostics now
operate on local slabs, as do analytic coarse/fine initialization and checkpoint read/write; see the [AMR case guide](../cases/amr_reactive_3d/README.md).

The regular 3D patch-tree geometry/state foundation now supports multiple children,
recursive PCM initialization, bottom-up restriction and composite inventories.
Snapshot peer/ancestor ghost filling and physical boundaries are connected.
Per-donor-level temporal interpolation is available for coarse ghosts.
Tree flux synchronization/reflux now handles nested and adjacent boxes, periodic
seams and physical boundary-ledger replacement. Synchronous SSPRK2 hydro now
connects flux production and temperature recovery, using a common timestep
bounded by all node CFL limits. Limited coarse-donor interpolation is connected
for PLM hydro and transport. Single-depth RK2 stepping now exposes integrated
fluxes and coarse-time endpoints. Interval-scoped reflux can return corrected
fluxes while preserving coarser endpoints. Recursive hydro subcycling and scaled
coarse CFL selection are connected. Recursive transport uses ratio-squared substeps
and scaled diffusion/wall bounds. Coupled/application dispatch and mode-aware
checkpoint continuation are now connected. MPI patch-owner scatter/gather is
implemented, and sparse owner-cell queries provide the payload exchange primitive;
synchronous halo donor selection, PCM/PLM interpolation and physical transforms
are connected. Distributed bottom-up restriction with depth limits is implemented.
MPI shared/periodic face flux synchronization is implemented. MPI application subcycling dispatch and mode-aware checkpoint identity are now connected.
Focused hydro-only scheduled/adaptive recursive execution and changed-rank restart
are checked. Optional deterministic cell-step-weighted whole-patch partitioning is
connected at initialization/regrid/restart. Final CSV now uses owner-local MPI-IO;
initial inventory and restart reference scales also avoid field gathers.
Checkpoint saves now use owner-local MPI-IO.
Checkpoint loading now streams an index with bounded root memory and reads node
ranges on their owners; optional indexed V2 eliminates the serial payload scan.
Direct same-geometry owner migration now transfers U/T without interpolation or
a root gather. A cost-driven rebalance API now accepts owner-reported patch costs
and migrates only for a sufficient predicted peak-load reduction. Scheduled recursive full-split execution now has a focused combined regrid,
rebalance and V2 changed-rank restart check. Temperature-triggered full-split
refinement also has this combined check. Particle/AMR application integration is
connected; large-tree scaling and broader EB/particle qualification remain open.
Owner-local kernel timing
and opt-in application rebalance cadence are now connected. Recursive MPI
reaction/transport/hydro composition is connected. MPI recursive transport now connects
ratio-squared steps, interval reflux and scaled diffusion/wall bounds. MPI recursive
hydro now connects ratio-based steps, interval reflux and scaled coarse CFL.
MPI transport level stepping with temporal donors, physical walls and integrated
fluxes is now connected. MPI hydro level stepping
with coarse-time donors and integrated face fluxes is connected. Sparse coarse-time ghost endpoints
and per-depth interpolation fractions are now connected. MPI temperature-adaptation application dispatch
and adaptation-aware changed-rank continuation are connected.
MPI temperature tagging now generates deterministic trees from owner-local hot
exposed cells; fixed/selected application controls and restart policy are connected.
Scheduled MPI tree remapping and event-aware checkpoint continuation are connected.
Sparse conservative MPI tree remapping now supports changed geometry/ratios/owners
with collective rollback and reactive temperature recovery.
Fixed/selected static MPI tree restart now supports changed rank counts and serial interchange.
Fixed/selected-mechanism static MPI tree application dispatch and final CSV output are connected. Distributed
transport/wall timestep selection and the full common-step operator composition
are connected. Synchronous molecular-transport RK2 is connected with moving/isothermal
physical walls and serial PCM/PLM donor parity. Exposed-cell
MPI chemistry and reaction/hydro Strang composition are connected. Synchronous
distributed SSPRK2 hydro is connected with serial PCM/PLM parity. Owner-local EOS recovery
and global minimum CFL selection are implemented. Parent flux
averaging, conservative cell correction, physical boundary accounting and depth-limited
distributed reflux are connected. Broader qualification remains. Exposed-cell chemistry and reaction/hydro splitting are connected;
synchronous molecular transport now also connects nested reflux and physical
walls. The transactional reaction/transport/hydro split API is connected.
Conservative overlap remapping now handles topology/refinement-ratio changes.
The serial application now connects static topology input, analytic initialization,
split evolution, scheduled conservative regridding and exposed-cell output.
Configure-time selected-mechanism dispatch now uses an isolated tree runtime.
Tree snapshot I/O, application stop/restart and schedule-position checks are implemented. Temperature-driven tiled topology generation, successful-step application cadence
and automatic-geometry restart are connected for serial and MPI ownership. See the
[tree design](design_decisions/amr-patch-tree-3d.md).

Particle Cartesian MPI migration preserves parcel/component payloads and local
ID ordering, with collective input/capacity rejection. Distributed gas feedback,
thermodynamic stepping, scheduled injection, breakup and rank-changing restart
are connected to applications. AMR film ownership, regridding/rebalancing and
restart are also connected
([dated film evidence](validation/spray-amr-film-application-20260918.md)). These
implementation results do not establish experimental spray/film validation,
production scaling or complete reference-solver parity.

## Contract for every increment

Before modifying code, specify the failure/capability, supported dimensions,
physics/parallel mode and a measurable exit test. Include a focused unit test,
an application/conservation or external-parity test, and a precise scope
statement. Separate numerical changes from structural refactoring.

For comparisons, record reference revision, mechanism/source hash, initial
state, observable, tolerance, comparison time and mesh/time-step convergence.
Do not loosen tolerances or regenerate golden outputs merely to obtain green
CI. Failed, interrupted, missing and completed checks are distinct states.
New work starts from the integrated main; a deliberate PR dependency must
state its base and merge order rather than recreating an accidental stack.

## Historical detail

The old milestone-by-milestone roadmap remains at the
[immutable checkpoint](https://github.com/GTT10/Fortranslate/blob/509d07929d52d9e570a6af1f6ec7e5836bd4fe5d/docs/completion_roadmap.md).
Detailed rationale remains in [design decisions](design_decisions/), the
[historical implementation inventory](implementation_status.md) and dated
[validation records](validation/). [Issue #101](https://github.com/GTT10/Fortranslate/issues/101)
tracks completion and actual tooling limitations. The unexecuted ancestor-only
branch-deletion audit is an administrative remainder, not a numerical feature.
