# Current status: PeleF 0.246.0

This page separates implemented capability, dated reported results, and
independent checks. The preceding integration is
recorded in [PR #102](https://github.com/GTT10/Fortranslate/pull/102). The new
opt-in spray scope and qualification boundary are in
[the spray guide](spray.md) and [0.246.0 validation](validation/0.246.0.md).
Older milestone prose remains historical detail, not a current global status.

## Supported boundary

Regular fixed/selected gas input now selects periodic SDC with
external forcing through `time_integrator='sdc'`
([checks](validation/sdc-app-20260926.md)). Molecular transport now participates
as an old/new source family ([checks](validation/sdc-transport-20260926.md)).
The library couples chemistry into periodic SDC
([checks](validation/sdc-reactions-20260926.md)), and fixed/selected application
input now enables it with persistent reaction sources
([checks](validation/sdc-reacting-app-20260926.md)). Broader boundary/AMR/MPI
integration remains open.
The SDC library additionally supports zero-gradient outflow with separate
boundary inventories ([checks](validation/sdc-outflow-20260926.md)). Fixed/selected
application routing now exposes this path, including a nonuniform reacting
check ([evidence](validation/sdc-outflow-app-20260926.md)). Outflow molecular
transport now participates with separate boundary accounting
([checks](validation/sdc-outflow-transport-20260926.md)); other physical
boundary types remain open.
Prescribed outlet pressure now reaches both SDC hydro and transport
([checks](validation/sdc-pressure-outlet-20260926.md)); existing ghost-state
admissibility limits still apply.
Specified-state and Mach inflow now reach SDC hydro and molecular transport
([checks](validation/sdc-inflow-20260926.md)). Slip, no-slip and isothermal
walls now connect reflection, viscous stresses and heat transport
([checks](validation/sdc-walls-20260926.md)).

The serial AMR library now runs common-step SDC across the whole patch tree,
including hydro, molecular transport, forcing, reactions, reflux and final
restriction ([checks](validation/sdc-tree-loop-20260926.md)). All patches finish
one stage before the next begins. Face-flux export and caller-owned halo APIs
remain available ([fluxes](validation/sdc-flux-export-20260926.md),
[halos](validation/sdc-halo-callbacks-20260926.md)). Static serial fixed/selected AMR applications now expose this path with separate
external-source inventories ([checks](validation/sdc-tree-app-20260926.md)).
Serial SDC checkpoint/restart now preserves reaction estimates and forcing
inventories in plain/indexed records ([checks](validation/sdc-tree-restart-20260926.md)).
Scheduled/automatic serial SDC regridding now transfers signed reaction estimates
with the gas ([checks](validation/sdc-tree-regrid-20260926.md)).
SDC hydro stages now support level selection and separate coarse temporal
donors for state and predictor sources ([checks](validation/sdc-level-hydro-20260926.md)).
Transport stages now select levels and interpolate parent endpoints; reaction
stages update one level while retaining other levels' state and source estimates
([checks](validation/sdc-level-sources-20260926.md)). Coupled level SDC now advances
one interval using parent state/source endpoints and exports integrated fluxes
([checks](validation/sdc-level-loop-20260926.md)). Serial recursive SDC now connects refinement-ratio substeps, parent temporal
donors, integrated-flux synchronization and application restart
([checks](validation/sdc-subcycle-20260926.md)). MPI now has an owner-local SDC
reaction stage with collective rollback and optional level selection
([checks](validation/sdc-mpi-reaction-20260926.md)). MPI SDC hydro now communicates
signed predictor sources and returns owner-local rates/integrated fluxes
([checks](validation/sdc-mpi-hydro-20260926.md)). MPI transport also returns
owner-local SDC rates/fluxes with collective limiter reduction
([checks](validation/sdc-mpi-transport-20260926.md)). Common-step MPI SDC now
couples hydro, transport, forcing, reactions and reflux while retaining owner-local
state ([checks](validation/sdc-mpi-loop-20260926.md)). The MPI checkpoint library
now reads/writes SDC V3/V4 records with owner-local source payloads and changed-rank
interchange ([checks](validation/sdc-mpi-checkpoint-20260926.md)). MPI fixed/selected applications now expose common-step SDC, source-preserving
regridding/rebalancing and changed-rank restart
([checks](validation/sdc-mpi-app-20260926.md)). MPI SDC hydro now selects a level
and samples separate parent state/source temporal endpoints
([checks](validation/sdc-mpi-level-hydro-20260926.md)). MPI transport now also
selects levels and interpolates parent endpoints
([checks](validation/sdc-mpi-level-transport-20260926.md)). The coupled MPI level
SDC driver now uses these stages and exports interval fluxes and predictors
([checks](validation/sdc-mpi-level-loop-20260926.md)). Recursive MPI SDC now connects child substeps, accumulated-flux synchronization
and application dispatch ([checks](validation/sdc-mpi-subcycle-20260927.md)). EB now has a covered-cell-aware forced SDC reaction stage
([checks](validation/sdc-eb-reaction-20260927.md)); EB transport also exposes a transactional source/flux stage
([checks](validation/sdc-eb-transport-20260927.md)). EB CTU prediction now accepts conserved SDC source rates
([checks](validation/sdc-eb-hydro-source-20260927.md)). EB hydro now exposes conservative rates, integrated fluxes and boundary/wall
ledgers ([checks](validation/sdc-eb-hydro-stage-20260927.md)). Physical boundary filling now feeds the EB SDC hydro stage
([checks](validation/sdc-eb-boundary-stage-20260927.md)). Small-cell redistribution now has a transactional source-rate interface
([checks](validation/sdc-eb-redist-20260927.md)). The serial single-grid EB SDC loop now couples hydro, transport, forcing,
reactions and redistribution ([checks](validation/sdc-eb-loop-20260927.md)).
Fixed/selected single-grid applications now accept SDC input and report
separate forcing inventories ([checks](validation/sdc-eb-app-20260927.md)).
Single-grid EB SDC checkpoint/restart now preserves reaction estimates and
forcing inventories ([checks](validation/sdc-eb-restart-20260927.md)).
The EB tree now exposes provisional SDC reaction updates across all levels or
a selected level, preserving refined-parent estimates until synchronization
([checks](validation/sdc-eb-tree-reaction-20260927.md)).
Signed EB tree SDC source halos now support peer/ancestor donors, ancestor
time interpolation and scalar physical extrapolation
([checks](validation/sdc-eb-tree-source-halo-20260927.md)).
The common-step EB tree now exposes source-aware provisional hydro rates and
time-integrated fluxes without averaging, redistribution or reflux
([checks](validation/sdc-eb-tree-hydro-20260927.md)).
Provisional EB tree transport now returns rates, integrated fluxes and a
limiter without averaging or synchronization
([checks](validation/sdc-eb-tree-transport-20260927.md)).
The serial EB tree now couples those stages in a common-step SDC loop, with
final peer/coarse-fine synchronization and separate source inventories
([checks](validation/sdc-eb-tree-loop-20260927.md)).
The serial EB tree checkpoint envelope now preserves SDC forcing inventories
and reaction fields, with transactional reads
([checks](validation/sdc-eb-tree-checkpoint-20260927.md)).
The serial tree application now routes static common-step SDC input and
checkpoint/restart, with an independent forcing budget column
([checks](validation/sdc-eb-tree-app-20260927.md)).
EB tree SDC regridding now transfers signed reaction estimates together with
gas state; scheduled coarsen/refine restart is checked
([checks](validation/sdc-eb-tree-regrid-20260927.md)).
EB SDC hydro and transport stages now support a selected level with ancestor
time interpolation; hydro separately interpolates source predictors
([checks](validation/sdc-eb-tree-level-stages-20260927.md)).
The coupled EB tree SDC loop now also advances a selected level, retains
inactive levels and exports interval fluxes plus old/new source predictors
([checks](validation/sdc-eb-tree-level-loop-20260927.md)).
Recursive EB tree SDC now accumulates interval fluxes, synchronizes adjacent
levels and routes through the serial application with checkpoint/restart
([checks](validation/sdc-eb-tree-subcycle-20260927.md)).
MPI EB trees now expose an owner-local SDC reaction stage, including selected
levels and collective failure handling
([checks](validation/sdc-mpi-eb-reaction-20260927.md)).
MPI EB transport now returns owner-local SDC rates and integrated fluxes while
preserving provisional parent states, with collective rollback
([checks](validation/sdc-mpi-eb-transport-20260927.md)).
MPI EB SDC hydro now exchanges signed predictor sources and returns owner-local
CTU rates and integrated fluxes, retaining provisional parent states
([checks](validation/sdc-mpi-eb-hydro-20260927.md)).
Common-step MPI EB SDC now couples reactions, old/new transport and external
sources, CTU hydro, final reflux and composite inventories
([checks](validation/sdc-mpi-eb-loop-20260927.md)). MPI EB SDC checkpoint V4 now stores owner-local reaction estimates and forcing
inventories with direct patch offsets, reads serial V3 records and supports
serial readers and changed-rank reads
([checks](validation/sdc-mpi-eb-checkpoint-20260927.md)). MPI EB application routing now supports fixed-topology common-step SDC with
rank-changing checkpoint continuation
([checks](validation/sdc-mpi-eb-app-20260927.md)). MPI EB SDC load rebalancing now migrates reaction estimates atomically with gas
and temperature ([checks](validation/sdc-mpi-eb-rebalance-20260927.md)). MPI EB SDC regridding now transfers signed reaction estimates with fluid-volume
restriction and PCM prolongation; scheduled changed-rank restart is checked
([evidence](validation/sdc-mpi-eb-regrid-20260927.md)). MPI EB SDC source halos now interpolate ancestor time endpoints with collective
control agreement and rollback
([checks](validation/sdc-mpi-eb-source-time-20260927.md)). MPI EB SDC hydro/transport stages now select a level and receive ancestor gas
endpoints; hydro additionally accepts separate source endpoints and fractions
([checks](validation/sdc-mpi-eb-level-stages-20260927.md)). The MPI EB SDC loop now advances one selected level, forwards distinct gas/source
times and exports fluxes plus old/new predictor fields
([checks](validation/sdc-mpi-eb-level-loop-20260927.md)). MPI EB SDC now recursively advances child levels, accumulates their fluxes and
synchronizes each parent interval. MPI application subcycling, coupled transport
timestep bounds and rank-changing continuation are connected
([checks](validation/sdc-mpi-eb-subcycle-20260927.md)). A serial composite flux-redistribution
operator transfers signed rates across peer, coarse/fine and paired periodic
leaf faces. Common-step serial SDC now applies it to hydro, molecular transport
and the final reflux increment
([operator](validation/sdc-eb-tree-redistribution-20260927.md),
[stage checks](validation/sdc-eb-tree-redistribution-loop-20260927.md)).
Serial application input and checkpoint identity now expose `flux_redist`
for common-step SDC ([checks](validation/sdc-eb-tree-redistribution-app-20260927.md)).
MPI common-step SDC and application dispatch now use owner-local composite
flux redistribution through sparse neighbor-rate and excess exchanges;
rank-changing continuation retains the selector
([checks](validation/sdc-mpi-eb-redistribution-20260927.md)). Recursive tree state redistribution and MPI EB spray SDC remain open. A serial
composite weighted-state operator now supports outward receivers and connected
neighborhood fallback across coarse/fine leaf faces, with conservation and
single-grid comparisons
([checks](validation/eb-tree-state-redistribution-20260927.md)); serial common-step SDC and application input now expose `state_redist` with
checkpoint-bound target fraction
([checks](validation/sdc-eb-tree-state-redistribution-20260927.md)). MPI now matches
serial common-step StateRedist through sparse value/mean exchanges and supports
rank-changing application restart
([checks](validation/sdc-mpi-eb-state-redistribution-20260927.md)). Recursive
StateRedist remains unconnected. Serial/MPI
source-level redistribution now separates immediate rates from time-integrated
cross-level transfers with a checked composite reconstruction
([checks](validation/sdc-eb-redistribution-level-transfer-20260927.md)); the
selected-level serial/MPI SDC loop now emits the accepted final-iteration
transfer alongside integrated fluxes
([checks](validation/sdc-eb-redistribution-level-loop-20260927.md)). The serial recursive
driver now accumulates transfers until their coarser source/recipient interval
closes, redistributes reflux increments and exposes application restart
([checks](validation/sdc-eb-redistribution-recursive-20260927.md)). MPI recursive
interval registers and application dispatch now follow the same owner-local
schedule, with exact serial/MPI state comparison and rank-changing continuation
([checks](validation/sdc-mpi-eb-redistribution-recursive-20260927.md)).

The variable-density forced chemistry RHS is now available as a library
building block ([checks](validation/forced-reactor-rhs-20260926.md)). An adaptive
implicit substep now integrates these rates
([checks](validation/forced-reactor-integrator-20260926.md)). Transactional
full-interval subcycling is also implemented
([checks](validation/forced-reactor-interval-20260926.md)) and feeds this library coupling.

Regular fixed/selected gas applications now accept force-per-volume input
`external_forcing(3)` with PPM/PeleC, time-centered work and separate source
inventories ([checks](validation/external-force-app-20260926.md)).

The full twelve-solve PPM CTU assembly now agrees with compiled pinned PeleC
bodies in eighteen small-grid cases, including independent mixture variation,
active shock flattening, artificial viscosity and inflow/outflow/slip/no-slip
boundaries ([comparison](validation/ppm-pinned-boundary-20260926.md)). This
includes normal tracing, Riemann solves, transverse correction, inflow PLM
replacement and final species normalization. Complete solver and physical
validation remain open.

Optional primitive sources now reach the periodic and supplied-halo PPM CTU
flux APIs, with 28 compiled-reference assembly cases
([evidence](validation/ppm-ctu-source-20260926.md)). Cell-local conserved-source
conversion now follows the pinned policy
([checks](validation/ppm-source-conversion-20260926.md)). Conserved source grids
now feed periodic/halo CTU directly, with 38 compiled-reference assembly cases
([checks](validation/ppm-conserved-source-20260926.md)); application source
generation and integration remain open.

The pressure path preserves independent and signed face species weights and
uses contact-speed passive upwinding
([projection](validation/ppm-mixture-projection-20260926.md),
[transverse](validation/ppm-signed-transverse-20260926.md),
[face flux](validation/ppm-passive-flux-20260926.md)). Physical wall/outflow
Riemann corrections reach regular serial/MPI AMR and MPI spray dispatch, while
explicit pressure outlets retain their prescribed donors
([checks](validation/ppm-physical-riemann-20260926.md)).

EB PLM normal prediction also preserves independently limited species weights
at open faces, with final CTU species-flux normalization retained
([regression](validation/eb-plm-species-20260926.md)).

Serial single-level EB KH eligibility now previews the complete boundary flight,
including periodic wraps and reflecting domain walls. Fixed-mechanism periodic
and wall cases preserve children and pending inventories through restart
([checks](validation/serial-eb-boundary-khrt-20260926.md)).

The serial EB impact interval now supports periodic crossings and specular
domain reflection before EB deposition, using the shared periodic film-image
counter ([checks](validation/serial-eb-periodic-impact-20260926.md)). Single-level
coupling and fixed/selected application restart now connect periodic film
exchange and reflecting-domain EB impacts
([evidence](validation/serial-eb-periodic-film-application-20260926.md)).

MPI EB deposition/splash and stationary films can now coexist with specular
outer-domain walls. Remaining flight and separate exterior/EB budgets are
preserved through regrid and changed-rank continuation
([checks](validation/mpi-eb-film-domain-reflection-20260926.md)). Periodic film thickness queries now sum all overlapping images across ranks
([checks](validation/mpi-periodic-film-query-20260926.md)). Periodic EB impact
flight and film/gas exchange now connect through regrid and changed-rank restart
([evidence](validation/mpi-eb-periodic-film-20260926.md)). Film creation on outer
domain faces remains open.

MPI EB KH child generation now accepts explicit nonperiodic domain walls,
with wall-before-release ordering and CTU regrid/changed-rank continuation
([checks](validation/mpi-eb-wall-khrt-20260926.md)). Explicit-boundary PCM
tree checkpoints now use a scoped embedded-record context, with changed-rank
continuation checked ([repair](validation/eb-tree-pcm-checkpoint-20260926.md)).

Patch Euler transport now accepts isothermal wall temperatures with transactional
heat-flux checks ([evidence](validation/eb-isothermal-patch-20260926.md)).
Serial common/recursive tree temperature forwarding now passes heated-wall
conservation checks ([evidence](validation/eb-tree-isothermal-transport-20260926.md)).
MPI common/recursive transport now agrees on wall temperatures and matches
serial ([checks](validation/mpi-eb-isothermal-transport-20260926.md)). Thermal
timestep bounds now connect in serial with an active-bound check
([evidence](validation/eb-tree-isothermal-timestep-20260926.md)); MPI thermal
selection now matches serial
([checks](validation/mpi-eb-isothermal-timestep-20260926.md)). Full reacting
common/recursive splits now forward wall temperatures with serial/MPI state
parity and rollback checks
([evidence](validation/eb-tree-isothermal-full-20260926.md)). Fixed/selected
serial and MPI gas applications now accept thermal wall inputs and bind wall
temperatures to changed-rank restart
([checks](validation/eb-tree-isothermal-application-20260926.md)). MPI EB
spray also forwards thermal gas walls, with common/recursive conservation and
fixed-mechanism changed-rank continuation checks
([evidence](validation/mpi-eb-isothermal-spray-20260926.md)).

MPI EB particle drift now supports specular domain-wall reflection with exterior
impulse accounting ([checks](validation/mpi-eb-wall-particle-20260926.md)).
Coupled common/recursive CTU wall steps and the fixed MPI spray application
now include wall reflection and changed-rank restart checks
([evidence](validation/mpi-eb-wall-spray-application-20260926.md)).

Fixed/selected tree gas applications now expose slip/no-slip adiabatic walls,
with scheduled regridding and changed-rank restart checked on an enclosed
planar fixture ([evidence](validation/eb-tree-wall-application-20260926.md)).

Serial EB tree halo construction now accepts adiabatic slip/no-slip walls,
including reflected six-layer corner donors with nonzero momentum
([checks](validation/eb-tree-wall-halo-20260926.md)). MPI wall halo exchange
now matches serial with distributed and empty-rank ownership
([checks](validation/mpi-eb-wall-halo-20260926.md)). Serial common PCM/CTU and
recursive PCM hydro now accept tangential adiabatic walls with composite
conservation checks ([evidence](validation/eb-tree-wall-hydro-20260926.md));
MPI common/recursive CTU wall hydro now matches serial for a tangential planar
fixture ([checks](validation/mpi-eb-wall-hydro-20260926.md)). Full reacting
wall splits now match serial for common/recursive CTU with explicit/implicit
chemistry ([checks](validation/mpi-eb-wall-full-20260926.md)). Serial common/recursive molecular transport now forwards
domain-wall labels and enforces impermeable adiabatic fluxes
([checks](validation/eb-tree-wall-transport-20260926.md)); MPI wall transport runtime is covered by those full-split checks.

Serial EB tree halos also accept paired periodic faces, wrapping donor indices
through peers and ancestors, including halos wider than the root domain
([checks](validation/eb-tree-periodic-halo-20260926.md)). MPI periodic tree
halo exchange now matches serial face-only, full one-layer and six-layer
construction, including temporal ancestor donors and idle ranks
([checks](validation/mpi-eb-periodic-halo-20260926.md)). Periodic tree application
dispatch is not yet connected.
The shared EB peer-flux correction kernel now accepts a periodic domain
length on one axis, checks matching face centroids and applies conservative
state correction across the seam
([checks](validation/eb-periodic-peer-flux-20260926.md)). Serial tree traversal
now reconciles periodic peers and mirrors restricted fine transfers into the
opposite coarse face. Common-step PCM and PLM/CTU hydro connect that path
([checks](validation/eb-tree-periodic-flux-20260926.md)). Serial common and
recursive molecular transport and recursive hydro now forward periodic
controls, including provisional endpoint correction before descendant time
interpolation ([checks](validation/eb-tree-periodic-recursive-20260926.md)).
Periodic common/recursive reacting PLM/CTU splits now pass focused explicit
and implicit chemistry checks. MPI flux synchronization also accepts periodic
faces with exact serial agreement, owner-local transfers and collective
rollback ([checks](validation/eb-tree-periodic-full-mpi-flux-20260926.md)).
Common and recursive MPI hydro now forward periodic labels through provisional
endpoint and final flux synchronization, with focused PLM/CTU serial agreement
and collective rollback checks
([checks](validation/mpi-eb-periodic-hydro-20260926.md)). MPI common/recursive
transport and reacting splits also connect periodic synchronization, with
serial agreement and collective rollback checks
([checks](validation/mpi-eb-periodic-full-20260926.md)). Fixed/selected serial
and MPI gas applications now accept paired periodic input, wrap geometry
padding and retain boundary identity across scheduled regridding and changed-rank
restart ([checks](validation/eb-tree-periodic-application-20260926.md)).
MPI EB spray also accepts periodic input, with fixed-mechanism 2/4-rank
comparison and 2-to-4-rank restart checks
([checks](validation/mpi-eb-periodic-application-20260926.md)). Periodic
EB impact/film combinations are also connected, as described above.
The composite EB particle stencil now accepts periodic axes and checks
visibility along segmented, wrapped paths to donor centroids
([checks](validation/eb-tree-periodic-stencil-20260926.md)). MPI particle field
sampling and gas deposition now route those periodic donors to their owners,
with conservation and rollback checks
([checks](validation/mpi-eb-periodic-field-20260926.md)). Mixed-particle drift
and reacting MPI EB spray coupling now forward periodic faces through field
exchange, drift and migration, with common/recursive CTU conservation checks
([checks](validation/mpi-eb-periodic-step-20260926.md)). Periodic KH child
generation now follows wrapped parent flight and passes changed-rank continuation
([evidence](validation/mpi-eb-periodic-khrt-20260926.md)).

Serial EB tree halo construction now accepts prescribed inflow states for
one-cell and six-cell halos while retaining interior peer/ancestor donors
([checks](validation/eb-tree-prescribed-halo-20260926.md)). Common and recursive
hydro/transport/reacting library updates now forward those conditions
([checks](validation/eb-tree-inflow-step-20260926.md)). MPI halo exchange also
accepts collectively agreed prescribed states, with serial agreement and idle-rank
rollback checks
([checks](validation/mpi-eb-prescribed-halo-20260926.md)). Common and recursive
MPI hydro now forward these controls with serial agreement
([checks](validation/mpi-eb-inflow-hydro-20260926.md)). MPI molecular transport,
full reacting splits and donor timestep bounds are also connected
([checks](validation/mpi-eb-inflow-full-20260926.md)). Fixed/selected serial and
MPI gas application inputs now connect prescribed inflow and boundary-bound
changed-rank restart across scheduled regridding
([checks](validation/eb-tree-inflow-application-20260926.md)). Tree application
inlets also accept per-face named mole fractions independently of initial gas,
with resolved composition bound to restart
([checks](validation/eb-tree-inlet-composition-20260926.md)). MPI EB spray now
also forwards prescribed inflow through the coupled split and application,
with fixed-mechanism changed-rank continuation checks
([checks](validation/mpi-eb-spray-inflow-20260926.md)).

Single-level fixed/selected EB gas applications now accept named mole fractions
per inlet face, with exact input-composition restart binding
([checks](validation/eb-inlet-composition-20260926.md)). Single-level spray
also constructs these inlet states and binds composition to its coupled restart
([checks](validation/eb-spray-inlet-composition-20260926.md)).

Single-level fixed/selected EB applications now offer `hydro_method='plm_ctu'`
with the `pelec` Riemann solver, outflow halos, reaction/transport composition
and method-bound checkpoint continuation. PCM remains the default. Fixed and
selected H2O2 coupled continuation passed focused checks
([evidence](validation/eb-ctu-application-20260925.md),
[input example](../cases/reactive_eb_3d/plm_ctu.nml)). This does not qualify a physical wall-flow solution. Explicit inflow/outflow input,
donor timestep bounds and boundary-bound restart are also connected for fixed and
selected single-level gas applications
([checks](validation/eb-inflow-application-20260925.md)). Fixed domain-wall inputs also
connect reflected geometry, thermal timestep limits and temperature-bound restart
([checks](validation/eb-wall-application-20260925.md)). Paired periodic faces also connect
wrapped geometry, coupled transport/CTU and changed-boundary restart rejection
([checks](validation/eb-periodic-application-20260925.md)).
MPI EB spray also forwards CTU through the gas/particle split, with common and
recursive library checks and fixed-mechanism changed-rank application continuation
([evidence](validation/mpi-eb-spray-ctu-20260925.md)). Serial EB spray fixed/selected
applications also connect CTU with method-bound continuation
([evidence](validation/serial-eb-spray-ctu-20260925.md)). The allocated serial
EB spray library also forwards explicit inflow/outflow gas boundaries and
escaped particle inventories
([checks](validation/eb-spray-open-coupling-20260925.md)). Fixed/selected serial
spray input, inlet timestep bounds and boundary-bound restart are connected
([application checks](validation/eb-spray-inflow-application-20260925.md)). Periodic
serial spray library coupling now wraps field exchange and parcel drift with
matching gas faces
([checks](validation/eb-spray-periodic-field-20260925.md)). Fixed/selected serial
periodic input and restart are also connected
([checks](validation/eb-spray-periodic-application-20260925.md)); periodic
impact/film dispatch now has the focused checks linked above.
Serial domain-wall gas conditions now connect specular particle reflection,
wall thermal timestep bounds and fixed/selected temperature-bound restart
([checks](validation/eb-spray-wall-20260926.md)).
Serial/MPI fixed and selected EB tree gas applications now also accept CTU,
including scheduled regridding and method-bound changed-rank continuation
([application checks](validation/eb-tree-ctu-application-20260925.md)).
Serial EB tree library operators support common and recursive CTU steps with
reaction/transport composition ([checks](validation/eb-tree-ctu-subcycle-20260925.md)).
MPI common and recursive tree hydro and coupled reaction/transport stepping
also support CTU with owner-local state and serial agreement
([common hydro](validation/mpi-eb-ctu-common-20260925.md),
[recursive/coupled checks](validation/mpi-eb-ctu-recursive-20260925.md)).

The pinned nozzle Mach-inlet state kernel is now available for six faces, with
interior composition and extrapolated or prescribed pressure
([state-kernel checks](validation/mach-inlet-state-20260925.md)). Serial hydro
now constructs its ghosts at both SSPRK2 stages
([hydro checks](validation/mach-inlet-hydro-20260925.md)). Molecular transport
also rebuilds Mach inlet states per stage
([transport checks](validation/mach-inlet-transport-20260925.md)). The full serial
reaction/transport/hydro composition forwards the same controls
([split checks](validation/mach-inlet-split-20260925.md)). Regular fixed/selected
input and per-step inlet timestep limits are connected
([input checks](validation/mach-inlet-input-20260925.md)). Regular Smagorinsky,
WALE, Vreman and Dynamic LES operators also accept the Mach controls
([LES checks](validation/mach-inlet-les-20260925.md)). Regular fixed/selected spray
now forwards Mach input through its coupled operators and binds restart context
([spray checks](validation/mach-inlet-spray-20260925.md)). Distributed wide Mach
halos and regular MPI Dynamic LES now match serial operators
([MPI LES checks](validation/mach-inlet-mpi-les-20260925.md)). The MPI spray
two-layer ghost entry also supports inviscid and viscous Mach donors
([ghost checks](validation/mach-inlet-mpi-spray-ghosts-20260925.md)). MPI hydro,
molecular transport, static LES and mixture splitting now forward the controls
([step checks](validation/mach-inlet-mpi-stepping-20260925.md)). Regular MPI spray
input, timestep limits and rank-changing restart are connected
([application checks](validation/mach-inlet-mpi-application-20260925.md)). The
abrupt-start timestep restriction from unused Dynamic filter padding is corrected
([startup checks](validation/dynamic-active-face-timestep-20260925.md)). The AMR
ghost sampler now reconstructs Mach inlet states after donor interpolation
([ghost checks](validation/mach-inlet-amr-ghosts-20260925.md)). MPI AMR ghost
exchange forwards the same controls with collective agreement and rollback
([MPI ghost checks](validation/mach-inlet-mpi-amr-ghosts-20260925.md)). Serial AMR
hydro now rebuilds Mach inlets for common and recursive level steps
([hydro checks](validation/mach-inlet-amr-hydro-20260925.md)). Distributed common
and recursive hydro now forward collectively agreed Mach controls
([MPI hydro checks](validation/mach-inlet-mpi-amr-hydro-20260925.md)). AMR
serial molecular transport also forwards Mach controls through common and
recursive stages ([transport checks](validation/mach-inlet-amr-transport-20260925.md)).
MPI molecular transport now also forwards the controls through common and
recursive stages ([MPI transport checks](validation/mach-inlet-mpi-amr-transport-20260925.md)).
Serial AMR LES selectors and common/recursive stepping also accept Mach controls
for all four models ([LES checks](validation/mach-inlet-amr-les-20260925.md)).
MPI AMR LES selectors and common/recursive stepping now forward collectively
agreed controls ([MPI LES checks](validation/mach-inlet-mpi-amr-les-20260925.md)).
Serial/MPI reacting transport/hydro composition now forwards Mach controls through
common and recursive steps ([split checks](validation/mach-inlet-amr-split-20260925.md)).
AMR spray coupling forwards the controls through the gas composition and both
LES half steps ([spray checks](validation/mach-inlet-amr-spray-20260925.md)).
Serial fixed/selected AMR gas application input, current-donor timestep limits
and restart checks are connected
([application checks](validation/mach-inlet-amr-application-20260925.md)). MPI gas
input and collective inlet bounds also support rank-changing and serial/MPI restart
([MPI application checks](validation/mach-inlet-mpi-amr-application-20260925.md)).
AMR spray input, timestep selection and context-bound restart are also connected
([spray application checks](validation/mach-inlet-amr-spray-application-20260925.md)).
This completes Mach dispatch for these regular and patch-tree gas/spray paths;
it does not establish physical nozzle validation or complete PeleC parity.

Periodic Dynamic LES is now connected to serial/MPI regular patch-tree stepping
and the MPI AMR reacting-spray application, including filter-bound restart.
Focused fixed-mechanism MPI parity and 2-to-4-rank continuation passed on
2026-09-25 ([runtime evidence](validation/mpi-dynamic-amr-runtime-20260925.md)).
Adaptive regridding/rebalancing with Dynamic LES and selected H2O2 restart now
also have focused execution evidence ([combined checks](validation/les-dynamic-adaptive-20260925.md)).
Prescribed reservoir inflow and extrapolated outflow also run through Dynamic
AMR LES, including changed-rank spray restart
([focused checks](validation/les-dynamic-open-20260925.md)).
Dynamic wall faces are connected through the AMR selector, stepper and spray application
([wall integration](validation/les-dynamic-amr-walls-20260925.md)).
Dynamic pressure-target outlets now correct all donor layers using physical
domain coordinates ([pressure integration](validation/les-dynamic-amr-pressure-20260925.md)).
Serial LES subcycling now has a recursive driver with interval-weighted reflux
([focused checks](validation/les-amr-subcycle-20260925.md)) and a
[scaled timestep selector](validation/les-subcycle-timestep-20260925.md).
MPI LES now also schedules recursive owner-local intervals and synchronizes
accumulated fluxes ([MPI recursion checks](validation/les-mpi-amr-subcycle-20260925.md)).
The MPI spray application now connects recursive LES half steps and the scaled
diffusion selector ([coupling checks](validation/les-mpi-spray-subcycle-20260925.md)).
The regular serial spray application also connects Dynamic physical-boundary
fluxes and restart ([wall/restart check](validation/regular-dynamic-spray-20260925.md)).
Regular MPI spray also connects Dynamic physical-boundary exchange, with
changed-rank wall/restart evidence
([MPI application check](validation/mpi-dynamic-spray-walls-20260925.md)).
Detailed-fuel Dynamic LES and physical qualification remain open.
Non-regular EB LES is an extension beyond the pinned reference, which rejects
such cells in all four closures ([source audit](validation/les-reference-scope-20260925.md)).
Use the regular-grid AMR LES path for the supported LES scope.

| Area | Implemented subset | Still open |
| --- | --- | --- |
| Reactive flow and transport | 1D/2D paths; regular 3D NASA7 hydro/chemistry/molecular transport with per-face periodic, prescribed reservoir inflow, outflow, slip-wall and adiabatic/isothermal no-slip boundaries with prescribed tangential wall motion | Characteristic inflow and moving geometry, general 3D boundary/coupling qualification |
| 3D AMR | Serial/MPI two-level and multi-patch, multi-level trees; prescribed or temperature-tagged regridding, conservative remap, physical boundaries, common or recursive reaction/transport/hydro, measured-work rebalancing and rank-neutral restart ([tree design](design_decisions/amr-patch-tree-3d.md)) | Broader EB/particle qualification, large-tree scaling and convergence |
| 3D embedded boundary | Planar, extruded-section and 3D nodal level-set single-level applications with fixed/selected reaction/transport/hydro, outward-neighbor StateRedist with connected-small-cell fallback and geometry-bound restart ([redistribution evidence](validation/eb-cluster-redistribution-20260917.md), [application evidence](validation/eb-level-set-application-20260917.md)) | Broader small-cell stabilization, convergence and integrated EB-AMR-MPI qualification |
| Chemistry | Fixed H2/O2 and configure-time selected supported bundles, including Lindemann/Troe/SRI/Tsang falloff, chemically activated rates, PLOG and Chebyshev; official-Fortran CVODE for selected 0D and opt-in regular 3D spray cells; FFCM1_Red methanol, GRI-Mech and Lu n-dodecane fixtures; input/runtime capacity up to 128 species | Diesel ignition/physical validation across conditions and sparse/performance qualification |
| Parallelism/output | Rank-parity and changed-rank restart for qualified subsets; distributed two-level MPI AMR CSV output | Production scaling, broader distributed I/O qualification and other output paths |
| LES and spray | Opt-in regular 3D single-component parcels with per-port liquid identity or shared or per-port multicomponent compositions, density/pressure fits and vapor-state caps, and periodic/reservoir-inflow/outflow/specular parcel boundaries and adiabatic/isothermal tangentially moving gas walls; temperature-dependent liquid density/heat capacity and per-port breakup viscosity/surface tension, selectable Schiller-Naumann or PelePhysics drag with ETAB deformation correction, drag/heating/evaporation with Clausius or bounded Antoine saturation pressure and an optional near-boiling vapor-state cap, and ETAB and KH-RT representative-drop breakup, optional KH shed-mass children with conservative transverse velocities, and bounded computational splitting, conservative two-way coupling and uniform external acceleration with impulse/work ledger, deterministic cone/scheduled injection from point or finite circular solid/hollow outlets with swirl with up to 16 independently configured ports and optional piecewise-linear shared or independent per-port injection-rate histories with optional discharge-based injection speed, and per-port mass- or number-CDF Weibull/normal/lognormal/ChiSquared size distributions, deviatoric Smagorinsky, coupled restart | PeleMP latent/boiling-property fidelity, physical mixture boiling, general dispersion/collision, broader characteristic-boundary and near-wall turbulence qualification, calibrated properties and experimental validation; lateral film transport is a future extension |

MPI AMR particle/gas feedback, film ownership and changed-rank application
restart are implemented ([film integration](validation/spray-amr-film-application-20260918.md));
MPI EB spray also connects transactional regridding and continuation
([EB integration](validation/spray-mpi-eb-regrid-application-20260918.md)). These
are historical focused results, not newly rerun qualification. A pinned-source
[property audit](validation/spray-property-scope-20260925.md) identifies the
remaining PeleMP latent/boiling distinction; GCM is restricted upstream to
PeleLMeX. The current constant-latent surface cap remains an approximation.

A library API now coarsens 3D EB geometry while preserving physical moments and
face closure; it rejects disconnected or unrepresentable coarse cells
([evidence](validation/eb-geometry-coarsening-20260917.md)). Reactive volume-weighted restriction and PCM prolongation now preserve conserved
components and recover EOS temperatures transactionally
([evidence](validation/eb-state-transfer-20260917.md)). Open-area flux restriction
and interior coarse/fine reflux now have conservative transactional library APIs
([evidence](validation/eb-reflux-20260917.md)). EB Riemann fluxes now accept face
halos, and plain Euler returns accepted flux integrals
([evidence](validation/eb-face-halo-20260917.md)). Geometry extraction and PCM
coarse-donor/outflow face halos now feed fine-patch Euler updates
([evidence](validation/eb-patch-halo-20260917.md)). A serial one-patch, two-level
Euler hydro API now connects both updates, reflux, averaging and separate composite
boundary/wall inventories transactionally
([evidence](validation/eb-two-level-hydro-20260917.md)). Masked explicit/implicit
reaction halves now wrap that update with composite chemical inventories
([evidence](validation/eb-two-level-chemistry-20260917.md)). Molecular-transport
fluxes now accept selectable face donors, including the species limiter
([evidence](validation/eb-transport-halo-20260917.md)). A synchronized two-level
transport step now connects limited fluxes, reflux and averaging
([evidence](validation/eb-two-level-transport-20260917.md)). A combined
reaction/transport/hydro step now accumulates all inventories transactionally
([evidence](validation/eb-two-level-full-20260917.md)). A fixed-mechanism serial
[application case](../cases/reactive_eb_two_level_3d/README.md) now runs this path
and writes composite state/inventory CSVs
([evidence](validation/eb-two-level-application-20260917.md)). Existing validated
EB checkpoint records can now share a caller-owned formatted stream
([evidence](validation/eb-embedded-checkpoint-records-20260917.md)). The two-level
application now saves both patches and cumulative budgets together, with staged
restore and exact focused continuation
([evidence](validation/eb-two-level-restart-20260917.md)). Rejected application
steps now retry with bounded halving and transactional state/budget preservation
([evidence](validation/eb-two-level-retry-20260917.md)). A plain Euler hydro
CFL library estimate now includes cut-cell fluid volume, open faces, embedded
wall area and optional remote signal bounds
([evidence](validation/eb-euler-cfl-20260917.md)). A plain Euler diffusion
interval now scales open-face conductance by fluid volume and accepts a remote
diffusivity bound ([evidence](validation/eb-patch-diffusion-cfl-20260917.md)).
The application now optionally selects a common interval from both grids,
with donor bounds and mode-aware restart
([evidence](validation/eb-two-level-auto-dt-20260917.md)). A library hydro path
now advances ratio-based fine substeps with coarse-time face donors, accumulated
fluxes and final reflux ([evidence](validation/eb-subcycle-hydro-20260917.md)).
The transport library now uses ratio-squared fine steps and accumulates
accepted limited fluxes ([evidence](validation/eb-subcycle-transport-20260917.md)).
Reaction halves now wrap both subcycled operators, and the application connects
scaled interval selection and mode-aware restart
([evidence](validation/eb-subcycle-full-20260917.md)). Same-level peer face
halos now have a transactional overlay API with shared-metric checks
([evidence](validation/eb-peer-halo-20260917.md)). An owner/receiver peer-flux
API now corrects receiver inventories and reconciles shared physical transfer
transactionally ([evidence](validation/eb-peer-flux-20260917.md)). A common-step
hydro library now advances multiple fine patches with deterministic face
ownership, peer synchronization and coarse reflux
([evidence](validation/eb-multipatch-hydro-20260917.md)). Molecular transport
now follows the same multiple-patch transaction with limited-flux reconciliation
([evidence](validation/eb-multipatch-transport-20260917.md)). Composite-masked
reaction halves now wrap the multiple-patch hydro/transport operators with
chemical inventories and a single state/diagnostic transaction
([evidence](validation/eb-multipatch-full-20260917.md)). The serial fixed-mechanism
application now connects static multiple patches, common-interval selection,
composite CSV output and transactional checkpoint continuation
([evidence](validation/eb-multipatch-application-20260917.md)). Multiple-patch
subcycling now synchronizes peer halos/fluxes on every fine step and accumulates
fluxes for final coarse reflux, with mode-bound restart
([evidence](validation/eb-multipatch-subcycle-20260917.md)); MPI EB remains open. A selected-mechanism executable now shares
the application path and binds bundle/composition identity on restart
([evidence](validation/eb-selected-two-level-20260917.md)). Generated bundles
now expose source element composition and a generic EB inventory API consumes
it ([evidence](validation/eb-element-metadata-20260917.md)). The two-level
selected frontend now uses generic elemental checkpoint records and chemical
conservation checks, including a 21-species methanol fixture
([evidence](validation/eb-carbon-checkpoint-20260917.md)). The single-level
selected frontend now shares source-derived element inventories and uses
schema-6 restart (fixed applications use schema 4). Redistributed hydro and
coupled APIs expose actual boundary, wall and chemical exchange vectors;
single-level applications accumulate and persist them, check all component
budgets and external-exchange-corrected elemental inventories. A density-contrast
methanol fixture passes exact state/budget restart
([application evidence](validation/eb-application-exchange-20260917.md),
[API evidence](validation/eb-redistributed-exchange-20260917.md),
[checkpoint evidence](validation/eb-exchange-checkpoint-20260917.md)).
A serial fixed-geometry two-level EB remapper now conservatively moves, adds
or removes fine patches, preserves same-resolution overlap, and supports
ratio changes with coarse-donor initialization and transactional failure.
Fixed/selected serial applications now connect prescribed patch schedules,
event-aligned stepping and schedule-bound checkpoint continuation; distributed
EB evolution remains open
([remap evidence](validation/eb-remap-20260917.md),
[scheduled application evidence](validation/eb-scheduled-application-20260917.md)).
An EB temperature-tagging API now selects deterministic coarse-index tiles
from exposed coarse/fine fluid temperatures, with buffered tagging and
transactional capacity rejection ([tagging evidence](validation/eb-temperature-tagging-20260917.md)).
Multipatch update APIs now accept zero fine patches and retain exact
single-grid split behavior and exchange budgets, enabling the numerical
endpoint for full derefinement ([coarse-only evidence](validation/eb-coarse-only-20260917.md)).
Fixed/selected applications now support coarse-only initialization, scheduled
2-to-0-to-2 patch transitions and exact checkpoint continuation through the
empty topology ([application evidence](validation/eb-zero-patch-application-20260917.md)).
Temperature-policy inputs now connect exposed-cell tagging, conservative
remapping and policy-bound adaptive restart at requested evaluation times.
High-temperature and no-hot-cell cases exercise both refinement and complete
derefinement ([adaptive evidence](validation/eb-adaptive-application-20260917.md)).
A core multilevel EB tree now reuses regular AMR topology for finest-lattice
geometry extraction, recursive PCM initialization, bottom-up restriction,
exposed-fluid inventories and transactional exposed-cell chemistry. Three-level
tests cover nested conservation and rollback; multilevel hydro/transport and
application integration remain open
([tree evidence](validation/eb-patch-tree-20260918.md)).
Spatial EB tree face halos now use a private averaged snapshot with same-level
peer precedence, ancestor fallback and physical outflow, including neighbors
under different parents. Multilevel flux synchronization and stepping remain
open ([halo evidence](validation/eb-tree-halo-20260918.md)).
A serial tree flux synchronization API now reconciles peer faces, restricts
fine transfers into both incident parent-level patches, corrects exposed
states and averages down in one transaction
([flux evidence](validation/eb-tree-flux-20260918.md)). Manufactured-transfer
conservation tests cover nested and one-sided refinement. A serial common-step
Euler hydro API now connects snapshot halos, all-node updates, synchronization
and separate exposed-wall/physical-boundary budgets transactionally
([hydro evidence](validation/eb-tree-hydro-20260918.md)). A matching molecular
transport API now uses limited node fluxes and globally identified closed
physical faces with transactional synchronization
([transport evidence](validation/eb-tree-transport-20260918.md)). A common-step
R/2 T/2 H T/2 R/2 API now combines these operators with separate boundary,
wall and chemical ledgers and complete rollback
([coupled evidence](validation/eb-tree-full-20260918.md)). A tree checkpoint
stream envelope now validates topology, all embedded node records and composite
budgets before transactional restoration; fixed-mechanism continuation is exact
in the focused fixture ([checkpoint evidence](validation/eb-tree-checkpoint-20260918.md)).
A serial fixed-mechanism
[tree application](../cases/reactive_eb_tree_3d/README.md) now connects static
input, coupled steps, exposed-cell/budget CSVs and staged checkpoint saving
with exact focused continuation ([application evidence](validation/eb-tree-application-20260918.md)).
Optional automatic common intervals now combine all-node cut-cell hydro and
diffusion bounds with global donor maxima and restart-bound policy
([interval evidence](validation/eb-tree-auto-interval-20260918.md)). MPI remains open for the EB tree application.
A conservative tree remap library now supports patch movement, coarsening and
refinement on an unchanged refinement lattice, with exact same-level overlap
and transactional geometry/EOS rejection
([remap evidence](validation/eb-tree-remap-20260918.md)). The application now
connects prescribed topology changes, exact event landing and schedule-bound
restart across coarsening/refinement ([schedule evidence](validation/eb-tree-scheduled-20260918.md));
a temperature-tagging library now reuses deterministic AMR tiling while
preserving exposed fine peaks and pruning solid-only nodes
([tagging evidence](validation/eb-tree-tagging-20260918.md)). The application now evaluates this policy at scheduled times, conservatively
remaps, and restores completed generated layouts under a policy-bound V4
checkpoint ([adaptive evidence](validation/eb-tree-adaptive-20260918.md)).
A selected-mechanism executable now shares the tree application, validates
composition and initial temperature range, checks chemical element conservation,
and binds bundle/composition identity on restart
([selected evidence](validation/eb-tree-selected-20260918.md)).
Tree face halos now optionally interpolate provisional ancestor endpoints
with independent per-level fractions while preserving current peer donors
([time-halo evidence](validation/eb-tree-time-halo-20260918.md)). A serial
recursive hydro API now advances ratio-based fine substeps, reconciles peers,
accumulates effective transfers and refluxes each parent interval
([subcycle hydro evidence](validation/eb-tree-subcycle-hydro-20260918.md));
a matching transport API now uses ratio-squared substeps and limited-flux
synchronization ([subcycle transport evidence](validation/eb-tree-subcycle-transport-20260918.md)).
The coupled API and application now select recursive hydro/transport inside
reaction splitting, with mode-bound scheduled restart
([subcycled application evidence](validation/eb-tree-subcycle-full-20260918.md)).
Automatic root intervals now scale patch hydro bounds by cumulative refinement
ratios and diffusion bounds by their squares
([interval evidence](validation/eb-tree-subcycle-timestep-20260918.md)).
An MPI library now restricts EB state/temperature on child owners and transfers
coarse payloads to parent owners, with replicated geometry validation and
exposed-cell composite reductions ([MPI transfer evidence](validation/mpi-eb-tree-transfer-20260918.md)).
MPI spatial face halos now fetch only requested owner-local U/T cells, using
same-level peers, ancestor PCM donors and physical outflow, with exact focused
serial parity ([halo evidence](validation/mpi-eb-tree-halo-20260918.md)).
MPI EB flux synchronization now reconciles peer faces, restricts and mirrors fine
transfers into parent peers, and refluxes exposed owner-local cells. A common
Euler hydro API connects these operations with face exchange and separate
boundary/wall ledgers, matching serial states exactly in focused cases
([MPI hydro evidence](validation/mpi-eb-tree-hydro-20260918.md)).
MPI common-interval molecular transport and exposed-cell reaction now wrap
hydro in transactional R/2 T/2 H T/2 R/2 composition, with focused serial state
parity for explicit/implicit chemistry ([MPI coupled evidence](validation/mpi-eb-tree-full-20260918.md)).
The fixed-mechanism `pelef_mpi_reactive_eb_tree_3d` application now connects
static or scheduled common-interval stepping, CSV/checkpoints and
changed-rank/serial restart ([MPI application evidence](validation/mpi-eb-tree-application-20260918.md)).
MPI automatic intervals now combine collective remote-donor maxima and owner-local
patch bounds, with exact focused serial timestep/state parity and changed-rank
restart ([automatic MPI evidence](validation/mpi-eb-tree-auto-20260918.md)).
Sparse MPI EB remapping now supports conservative patch movement, coarsening,
refinement and changed owners on the shared immutable lattice. Prescribed
application regrids and restart across their events are connected
([scheduled MPI evidence](validation/mpi-eb-tree-remap-20260918.md)).
MPI temperature adaptation now gathers exposed fluid-cell votes, preserves fine
hot spots and removes solid-only tiles. Adaptive application events and changed-rank
restart match the focused serial case
([adaptive MPI evidence](validation/mpi-eb-tree-adaptive-20260918.md)).
MPI ancestor halo donors now support independent time fractions at each depth,
with exact focused serial parity and collective rollback
([temporal MPI evidence](validation/mpi-eb-tree-time-halo-20260918.md)).
MPI recursive hydro now uses ratio-based fine intervals, provisional peer-corrected
endpoints and accumulated effective fluxes, matching focused serial U/T exactly
([MPI subcycle hydro evidence](validation/mpi-eb-tree-subcycle-hydro-20260918.md)).
MPI molecular transport now uses ratio-squared substeps and limited-flux
accumulation; the coupled API selects recursive H/T inside reaction splitting
([MPI recursive coupling evidence](validation/mpi-eb-tree-subcycle-full-20260918.md)).
The MPI application now dispatches recursive steps and mode-bound restart.
Automatic root intervals scale hydro bounds by cumulative refinement and diffusion
bounds by its square, with focused scheduled and automatic restart parity
([MPI subcycled application evidence](validation/mpi-eb-tree-subcycle-application-20260918.md)).
The selected-mechanism MPI EB executable now shares composition/element checks,
stepping controls and identity-bound checkpoints with serial. Its dedicated
runtime excludes fixed mechanism libraries, with focused H2/O2 implicit
changed-rank/serial restart parity
([selected MPI evidence](validation/mpi-eb-tree-selected-20260918.md)).
EB tree state CSV now uses owner-local MPI-IO chunks in serial row order,
with fixed/selected byte parity and preflight output protection
([distributed CSV evidence](validation/mpi-eb-tree-csv-20260918.md)).
Checkpoint saving now stages serial patch records on their owners and writes
them directly with MPI-IO; fixed/selected patch bytes and changed-rank/serial
restart match focused references
([distributed save evidence](validation/mpi-eb-tree-save-20260918.md)).
Checkpoint loading now scans record offsets on rank zero and restores each patch
only on its owner, with collective state/ledger publication and focused
changed-rank/serial parity
([distributed load evidence](validation/mpi-eb-tree-load-20260918.md)).
Initial tree fields now allocate only on patch owners, with sparse parent-section
fetches and the serial PCM kernel, removing full-tree state replication and scatter
([owner-local initialization evidence](validation/mpi-eb-tree-initialize-20260918.md)).
Root seed U/T is now constructed only on the root patch owner and released on
return from initialization, including restart topology setup
([evidence](validation/mpi-eb-tree-seed-20260918.md)). MPI checkpoints now include
a V2 patch byte-length index, allowing direct owner reads without a full root
scan. Serial readers also accept V2; legacy V1 saves retain the bounded-memory
scan fallback ([evidence](validation/mpi-eb-tree-index-20260918.md)). Geometry
remains replicated and scaling qualification remains open.
The MPI EB application can now rebalance whole patches using sampled owner-local
reaction/hydro/transport kernel time, with a configurable interval and minimum
estimated load improvement. Accepted common or recursive steps provide samples;
failed retry samples are discarded. State and temperature migrate together
without interpolation or a root gather
([evidence](validation/mpi-eb-tree-rebalance-20260918.md)). Timing-dependent owner
maps are performance state and do not constrain checkpoint restart identity.

A new regular 3D patch-tree library represents multiple nonoverlapping children
and arbitrary depth, with parent-before-child indexing and a common refinement
ratio at each depth. It reuses existing PCM/restriction kernels and provides
transactional initialization, bottom-up averaging and composite inventories
([design](design_decisions/amr-patch-tree-3d.md),
[evidence](validation/amr-patch-tree-3d-20260917.md)). This is a geometry/state
foundation. Snapshot ghost filling now exchanges peer data across parent boxes,
uses PCM ancestor fallback and follows regular-grid periodic/physical boundary
conventions ([evidence](validation/amr-tree-ghost-3d-20260917.md)). Per-level coarse
endpoint interpolation is also available for ghost donors while retaining exact
current peers ([evidence](validation/amr-tree-temporal-ghost-20260917.md)). Limited linear coarse-donor interpolation is now available for hydro
([evidence](validation/amr-tree-plm-ghost-20260917.md)) and transport
([evidence](validation/amr-tree-transport-plm-20260917.md)); coupled/application subcycle dispatch is opt-in. A single-depth RK2 hydro API now accepts coarse
time endpoints and returns integrated face fluxes without changing inactive levels
([evidence](validation/amr-tree-hydro-level-20260917.md)). Reflux now supports
completed intervals restricted to a minimum depth and returns effective fluxes
for parent accumulation ([evidence](validation/amr-tree-interval-reflux-20260917.md)). A recursive
hydro API now connects ratio-based substeps, coarse-time donors, interval reflux,
temperature recovery and a scaled coarse CFL bound
([evidence](validation/amr-tree-subcycle-hydro-20260917.md)). Recursive transport now
uses ratio-squared substeps, coarse-time donors, interval reflux and scaled
diffusive/isothermal-wall bounds
([evidence](validation/amr-tree-subcycle-transport-20260917.md)). Fixed/selected applications
now select the recursive reaction/transport/hydro composition, scaled timestep
bounds and mode-aware checkpoint continuation
([evidence](validation/amr-tree-subcycle-application-20260917.md)). MPI tree transfer APIs now
store payloads only on assigned patch owners, with collective metadata checks and
root scatter/gather ([evidence](validation/amr-tree-mpi-distribution-20260917.md)).
Same-geometry ownership changes now move conserved state and temperature directly
between owners with exact values and collective preflight rollback
([evidence](validation/amr-tree-mpi-redistribution-20260917.md)).
A cost-driven rebalance API combines owner-reported patch costs, a minimum
predicted gain and transactional migration. Fixed/selected MPI applications now
optionally measure owner-local kernels and rebalance at a configured step cadence
([evidence](validation/amr-tree-mpi-runtime-balance-20260917.md)).
Sparse cell queries now fetch only requested owner payloads in caller order, including
empty ranks and collective rollback ([evidence](validation/amr-tree-mpi-cells-20260917.md)).
Synchronous MPI halos now share serial peer/ancestor lookup, PCM/limited PLM and
physical boundary transforms ([evidence](validation/amr-tree-mpi-ghost-20260917.md)).
Final MPI tree CSV output now writes owner-local chunks directly with MPI-IO and
uses serial-order small-vector final inventory reduction
([evidence](validation/amr-tree-mpi-csv-20260917.md)). Initial inventory and the
restart initial absolute-inventory scale now use the same distributed traversal
([evidence](validation/amr-tree-mpi-initial-inventory-20260917.md)). Checkpoint
saves now validate EOS/parent averages collectively and write owner-local payloads
with serial-byte compatibility
([evidence](validation/amr-tree-mpi-checkpoint-write-20260917.md)). Restart now
reads indexed node ranges on their assigned owners and validates candidates
without constructing or scattering full root fields
([evidence](validation/amr-tree-mpi-checkpoint-read-20260917.md)). Optional V2
checkpoints now store validated relative offsets, eliminating that reader's
root payload scan while retaining V1 input support
([evidence](validation/amr-tree-indexed-checkpoint-20260917.md)).
MPI applications now optionally assign whole patches by estimated cell-step work
at initialization, regridding and restart
([evidence](validation/amr-tree-mpi-partition-20260917.md)).
MPI bottom-up restriction now sends restricted child boxes directly to parent
owners, with optional depth limits ([evidence](validation/amr-tree-mpi-restriction-20260917.md)).
MPI shared-face flux synchronization now exchanges overlapping face rectangles
and prioritizes fine-replaced values ([evidence](validation/amr-tree-mpi-flux-peers-20260917.md)).
Distributed reflux now connects parent face averaging, uncovered-cell correction,
physical boundary accounting and bottom-up state restriction, including depth-limited
intervals ([evidence](validation/amr-tree-mpi-reflux-20260917.md)).
Owner-local temperature recovery and global minimum CFL selection now reuse the
serial thermodynamic kernels ([evidence](validation/amr-tree-mpi-thermo-20260917.md)).
Synchronous MPI SSPRK2 hydro now connects these operations for owned patches,
with PCM/PLM and physical inflow/walls ([evidence](validation/amr-tree-mpi-hydro-20260917.md)).
Owner-local exposed-cell chemistry and reaction/hydro Strang splitting now support
explicit/implicit integration ([evidence](validation/amr-tree-mpi-chemistry-20260917.md)).
Synchronous MPI molecular-transport RK2 now connects viscosity, conduction, species
and pressure diffusion with physical/moving isothermal walls
([evidence](validation/amr-tree-mpi-transport-20260917.md)).
MPI transport/wall timestep selection and the full common-step reaction/transport/hydro
composition are now connected ([evidence](validation/amr-tree-mpi-split-20260917.md)).
The fixed-mechanism MPI tree executable now connects static topology, owner-local
initialization/stepping, timestep retries and serial-format final CSV/budget output
([evidence](validation/amr-tree-mpi-application-20260917.md)). MPI recursive hydro subcycling and its scaled coarse CFL bound are now connected
([evidence](validation/amr-tree-mpi-subcycle-hydro-20260917.md)). MPI recursive transport now connects ratio-squared steps, interval reflux and scaled
diffusion/wall bounds ([evidence](validation/amr-tree-mpi-subcycle-transport-20260917.md)).
MPI reaction/transport/hydro composition now accepts recursive stepping
([evidence](validation/amr-tree-mpi-subcycle-split-20260917.md)); fixed/selected MPI application mode dispatch, scaled timestep selection and mode-aware
checkpoint continuation are now connected
([evidence](validation/amr-tree-mpi-subcycle-application-20260917.md)).
Focused hydro-only scheduled/adaptive recursive cases now cover fixed/selected
serial parity and changed-rank restart
([evidence](validation/amr-tree-mpi-subcycle-dynamic-20260917.md)). Scheduled recursive full-split
execution now also combines chemistry, molecular transport, moving walls, PLM,
measured-work rebalancing and V2 changed-rank continuation with exact serial bytes
([evidence](validation/amr-tree-mpi-dynamic-split-20260917.md)). Temperature-triggered
full-split refinement now has the same combined V2/rebalance/changed-rank check,
including no repeated tag on restart
([evidence](validation/amr-tree-mpi-adaptive-split-20260917.md)). MPI transport now advances one depth with
coarse-time donors, moving/isothermal walls and integrated face fluxes
([evidence](validation/amr-tree-mpi-transport-level-20260917.md)). MPI hydro
now advances one depth with coarse-time donors and returns integrated face fluxes
([evidence](validation/amr-tree-mpi-hydro-level-20260917.md)). MPI ghost filling now
accepts coarse time endpoints and per-depth interpolation fractions, with optional
target-depth selection ([evidence](validation/amr-tree-mpi-temporal-ghost-20260917.md)).
Fixed/selected MPI applications now connect temperature-adaptation controls,
post-step tagging/remapping and adaptation-aware changed-rank checkpoint continuation
([evidence](validation/amr-tree-mpi-adaptation-20260917.md)).
A distributed temperature-tagging API now shares only hot exposed-cell boxes and
uses the serial tile builder to produce identical geometry across ranks
([evidence](validation/amr-tree-mpi-tagging-20260917.md)).
Fixed/selected MPI tree applications now connect scheduled sparse remapping, event-time
clipping, conservation checks and schedule-aware changed-rank restart
([evidence](validation/amr-tree-mpi-schedule-20260917.md)).
A sparse MPI conservative-remap API now transfers intersecting old exposed cells
to new owners, preserves serial accumulation order, and recovers temperatures
transactionally ([evidence](validation/amr-tree-mpi-remap-20260917.md)).
Fixed/selected MPI tree checkpoint/restart now reuses the serial format and validates
context/geometry/EOS/ledger before owner redistribution
([evidence](validation/amr-tree-mpi-restart-20260917.md)). The selected-mechanism
MPI tree frontend now uses an isolated selected runtime and forwards the bundle
integrator ([evidence](validation/amr-tree-selected-mpi-20260917.md)). A transactional tree
flux-correction API now reconciles shared and periodic faces, propagates nested
fine fluxes, corrects uncovered cells and returns the physical boundary-ledger
replacement ([evidence](validation/amr-tree-reflux-3d-20260917.md)). Synchronous
SSPRK2 hydro now connects PCM/PLM face fluxes, temperature recovery and the
composite boundary ledger, with a CFL minimum over all nodes
([evidence](validation/amr-tree-hydro-20260917.md)); PLM hydro now enables limited coarse donors.
Cell chemistry now runs on exposed tree cells, with bottom-up restriction and
transactional reaction/hydro Strang splitting
([evidence](validation/amr-tree-chemistry-20260917.md)). Synchronous molecular
transport now includes viscous/thermal/species fluxes, nested reflux and physical
wall budgets ([evidence](validation/amr-tree-transport-20260917.md)). A combined
R/2-T/2-H-T/2-R/2 update now aggregates those boundary budgets transactionally
([evidence](validation/amr-tree-split-20260917.md)). A conservative overlap-remap
API now supports changed patch layouts and refinement ratios, with transactional
EOS recovery ([evidence](validation/amr-tree-remap-20260917.md)). Automatic solution-based tagging and MPI ownership are now connected as described above. The experimental
[serial tree application](../cases/amr_tree_reactive_3d/README.md) now reads explicit
patch topology, initializes analytic fields, evolves the split system and writes
exposed-cell CSV plus composite budgets
([evidence](validation/amr-tree-application-20260917.md)). Prescribed-time topology
changes now connect conservative remapping, exact event-time stepping and budget
continuity ([evidence](validation/amr-tree-regrid-20260917.md)). Configure-time
selected mechanisms now use the same tree driver through an isolated runtime
([evidence](validation/amr-tree-selected-20260917.md)). Versioned tree snapshot I/O
and physics-context comparison APIs are now available
([evidence](validation/amr-tree-checkpoint-api-20260917.md)). Fixed/selected application
stop/restart now validates schedule position, geometry and continuation context
([evidence](validation/amr-tree-restart-20260917.md)). A temperature-tagged topology
builder now creates multiple tiled children over requested depths, including
periodic buffers and fine-only hot-cell preservation
([evidence](validation/amr-tree-tagging-20260917.md)). Successful-step adaptation
cadence and fixed/selected automatic-geometry restart are now connected
([evidence](validation/amr-tree-adaptation-20260917.md)).

Sparse MPI AMR has physical ghost layers and two-level PCM/PLM hydro supporting
one coarse column per rank, with serial equality and collective rollback checks
([ghosts](validation/amr-mpi-boundary-halo-20260917.md),
[coarse hydro](validation/amr-mpi-boundary-hydro-20260917.md),
[hierarchy](validation/amr-mpi-physical-hierarchy-20260917.md)). The
[reaction/hydro split](validation/amr-mpi-physical-reaction-20260917.md) also forwards
physical boundary data and budgets. Physical molecular transport and the full
reaction/transport/hydro split are also connected
([evidence](validation/amr-mpi-physical-transport-20260917.md)). MPI application
inputs, wall/inlet timestep limits and rank-neutral boundary-ledger restart are
also connected ([evidence](validation/amr-mpi-physical-app-20260917.md)). Fixed elementary
and full-H2/O2 reacting AMR checkpoint/restart now retains canonical reaction,
composition and integrator context for serial/MPI applications
([evidence](validation/amr-fixed-chemistry-restart-20260917.md)).
The selected AMR runtimes remain independent of fixed-mechanism databases and now
include distributed checkpoint I/O in their build
([build isolation](validation/selected-amr-checkpoint-isolation-20260917.md)).

Boundary-touching AMR now has a physical flux-correction API that replaces
covered domain-surface fluxes in the boundary ledger and corrects only uncovered
interior coarse cells ([evidence](validation/amr-boundary-reflux-20260917.md)).
Serial PCM/PLM hydro now connects fine physical ghosts, subcycling, boundary
ledger replacement and average-down for such patches
([evidence](validation/amr-boundary-touching-hydro-20260917.md)). Cell-local reaction,
reaction/hydro splitting and the state-based CFL helper now accept these patches
([evidence](validation/amr-boundary-touching-reaction-20260917.md)). Molecular
transport and the full reaction/transport/hydro split are connected, including
fine wall heat/work fluxes ([evidence](validation/amr-boundary-touching-transport-20260917.md)).
Serial static application inputs, fine wall timestep bounds and checkpoint/restart
are connected ([evidence](validation/amr-boundary-touching-app-20260917.md)).
Serial temperature-driven regridding now supports boundary-touching patches with
explicit physical tagging and geometry-aware restart
([evidence](validation/amr-boundary-regrid-20260917.md)). Sparse MPI distribution,
scatter/gather and cell-local reaction accept boundary-touching geometry
([evidence](validation/amr-mpi-boundary-distribution-20260917.md)). MPI PCM/PLM
hydro and reaction/hydro splitting now connect fine physical ghosts and covered
surface flux replacement ([evidence](validation/amr-mpi-boundary-touching-hydro-20260917.md)).
MPI molecular transport and the full reaction/transport/hydro split now also
connect fine wall heat/work and physical flux replacement
([evidence](validation/amr-mpi-boundary-touching-transport-20260917.md)).
MPI static application inputs now include fine wall timestep limits and
rank-neutral restart for boundary-touching geometry
([evidence](validation/amr-mpi-boundary-touching-app-20260917.md)).
MPI physical-boundary temperature tagging, conservative remap and changed-rank
restart are also connected
([evidence](validation/amr-mpi-boundary-regrid-20260917.md)). Temperature tagging
now scans local slabs and reduces only tag bounds; unchanged patches skip full
hierarchy gathering ([evidence](validation/amr-mpi-distributed-tagging-20260917.md)).
Conservative remapping now runs within fixed coarse slabs using neighboring
coarse halos and local fine overlap reuse
([evidence](validation/amr-mpi-distributed-remap-20260917.md)). It no longer
gathers full fields during regridding. An explicit repartition API now transfers
coarse/fine states and temperatures directly between ranks, and remapping
preserves the chosen coarse partition
([evidence](validation/amr-mpi-repartition-20260917.md)). Optional cell-count
balancing is connected after initialization/restart and geometry changes
([evidence](validation/amr-mpi-load-balance-20260917.md)). General patches/levels,
measured performance/scaling and broader physical qualification remain open. MPI AMR
CSV output now uses direct collective file writes, and final inventory, synchronization
and composition diagnostics use distributed reductions
([evidence](validation/amr-mpi-distributed-csv-20260917.md)). Analytic initial
conditions are also generated on local coarse/fine slabs, with local average-down
and distributed initial inventories
([evidence](validation/amr-mpi-local-initialization-20260917.md)). Checkpoint writes
now serialize the existing format directly from local slabs after collective
metadata, temperature and synchronization checks
([evidence](validation/amr-mpi-checkpoint-output-20260917.md)). Restart now reads
canonical checkpoint cell rows directly into local candidate slabs, then commits
geometry, fields and diagnostics only after collective checks
([evidence](validation/amr-mpi-checkpoint-input-20260917.md)). Metadata remains
replicated; reformatted cell records use a bounded-memory streaming fallback.
Production I/O scaling is not yet measured.

Ahamed wall-impact classification, splash velocities and azimuthal mass
fractions, four-child construction and deposited cylindrical film patches are
connected to the regular coupled-step API with wall budgets, stationary film
heating/evaporation with per-material constant or cubic thermal conductivity,
CLI inputs/output and coupled checkpoint persistence
[evidence](validation/wall-film-conductivity-20260916.md). Lateral film transport
and physical impact qualification remain open. The gas-side film skin now uses
the pinned half-temperature/half-composition rule
([correction evidence](validation/wall-film-skin-20260916.md)).

KH shed-mass accumulation and four-child release now connect to regular CFD
and restart ([evidence](validation/khrt-children-20260916.md)); computational
splitting now partitions the KH inventories
([evidence](validation/khrt-splitting-20260916.md)). Stationary wall impact
now consumes crossed parcels before KH release, with coupled film/lineage restart
([evidence](validation/khrt-wall-20260916.md)); wet-impact correlation limits remain.

The mixed-liquid API now covers mixture properties, component transport,
regular-grid feedback, component-specific boundary budgets, representative
breakup/splitting and composition checkpoint persistence
([evidence](validation/spray-mixture-grid-20260916.md)). The combined
hydro/transport/LES/chemistry driver now retains composition and rolls it back
with all other state ([evidence](validation/spray-mixture-reactive-20260916.md)).
Mixture CLI inputs, shared/independent scheduled injection, component CSVs and
restart compatibility are connected
([evidence](validation/spray-mixture-cli-20260916.md)). Trace-species depletion and
minimum-diameter removal now return their complete inventories to gas
([evidence](validation/spray-mixture-terminal-20260916.md)). Mixed KH children now inherit composition with conservative transverse velocities
and retain their state through splitting and restart
([evidence](validation/spray-mixture-kh-20260916.md)). The mixture impact event API now partitions component masses into splash and film
with wall budgets ([evidence](validation/spray-mixture-wall-event-20260916.md));
its CLI connection is now available. Planar crossing now uses mixture properties
and the reference mass-weighted impact boiling temperature
([evidence](validation/spray-mixture-wall-crossing-20260916.md)). A supplied film-exchange API now
updates composition, temperature and height with gas/wall energy closure
([evidence](validation/spray-mixture-film-exchange-20260916.md)). Mixed film
heat/evaporation rates and a frozen-rate update are now available
([evidence](validation/spray-mixture-film-rates-20260916.md)). The checkpoint API
now retains film component masses with transactional validation
([evidence](validation/spray-mixture-film-checkpoint-20260916.md)); application wiring is connected. A transactional mixed-film grid feedback API now
conserves species and gas/film/wall energy and removes terminal patches
([evidence](validation/spray-mixture-film-grid-20260916.md)). The mixed particle
coupler now connects wall crossing, splash/deposition and remaining-interval film
evolution in one transaction ([evidence](validation/spray-mixture-wall-coupled-20260916.md)).
The reactive driver now retains mixed film state across both spray half-steps and
rolls it back after late failure ([evidence](validation/spray-mixture-wall-reactive-20260916.md)).
CLI wall inputs, film component output, conservation totals and restart are connected
([evidence](validation/spray-mixture-wall-cli-20260916.md)). Active trace-film
components now transfer their full residual inventories to gas
([evidence](validation/spray-mixture-film-trace-20260916.md)). Per-port compositions now drive initial and
scheduled injection with component budgets and restart identity
([evidence](validation/spray-port-mixture-20260916.md)). Component density/heat-capacity fits, Antoine
pressure and Clausius surface caps are available in the CLI and restart context
([evidence](validation/spray-mixture-properties-20260916.md)).

Conservation, rollback, byte-parity, and short-time H2/O2 tests do not establish
physical ignition/spray validation or complete PeleC equivalence. Choosing a
license remains an owner decision.

AMR particle ownership now locates each parcel in the finest covering regular
3D tree cell and constructs stable owner-rank bins while preserving all parcel
fields and optional component inventories. Periodic lookup, half-open patch
faces, regrid reassignment and transactional validation are supported
([evidence](validation/spray-amr-partition-20260918.md)). MPI migration now uses
this ownership map and a shared Cartesian/AMR particle exchange, preserving
component and breakup payloads with collective validation/rollback through
particle crossings, owner changes and regridding
([evidence](validation/spray-amr-migration-20260918.md)). Serial composite-field
sampling/deposition now uses a shared positive cloud-overlap stencil over exposed
cells, with volume-scaled deposition and covered-parent averaging
([evidence](validation/spray-amr-coupling-20260918.md)). A serial gas-feedback
wrapper now recovers EOS temperatures after deposition and parent averaging,
committing gas U/T only when every cell accepts the update
([evidence](validation/spray-amr-gas-20260918.md)). MPI particle-position sampling
now requests only the composite stencil cells from their owners and matches
serial accumulation order ([evidence](validation/spray-amr-sampling-20260918.md)).
MPI field deposition now routes volume-scaled particle contributions to cell
owners and refreshes covered parents without a full-field gather
([evidence](validation/spray-amr-deposition-20260918.md)). MPI gas feedback now
recovers owner-local temperatures after parent averaging and commits U/T
collectively, with initial-state validation and failure rollback
([evidence](validation/spray-amr-mpi-gas-20260918.md)). A periodic MPI AMR
mixed-particle interval now connects drag/heating/evaporation, representative
breakup and KH children, common substeps, transactional gas feedback and final
finest-owner migration ([evidence](validation/spray-amr-step-20260918.md)).
An MPI reacting-spray step now composes particle half-intervals around the
existing reaction/transport/hydro split, with common or recursive gas steps and
whole-step rollback ([evidence](validation/spray-amr-reactive-20260918.md)).
A coupled remap API now commits tree/owners, conservative gas remapping and
finest-owner particle migration together, preserving all original state on
failure ([evidence](validation/spray-amr-remap-20260918.md)).
Fixed and selected [MPI AMR spray applications](../cases/mpi_amr_spray_3d/README.md)
now connect owner-local gas setup, mixed-cone initial injection, coupled stepping,
one scheduled regrid and gas/particle CSV output
([evidence](validation/spray-amr-application-20260918.md)). Explicit single-cone
injection schedules now land on event times and commit particle payloads, IDs and
cumulative injection inventories after successful migration
([evidence](validation/spray-amr-events-20260918.md)). AMR spray checkpoints now
bind gas/particle state, topology, IDs and event position, with immutable
publication and changed-rank continuation in fixed/selected executables
([evidence](validation/spray-amr-checkpoint-20260918.md)). Up to 16 AMR spray
ports now support per-port geometry, mixture fractions, nozzle profile/swirl and
equal-time injection events with preserved port IDs across restart
([evidence](validation/spray-amr-ports-20260918.md)). Temperature-driven AMR
adaptation now updates gas/particle topology on completed-step intervals and
continues after changed-rank restart without repeating the saved step's check
([evidence](validation/spray-amr-adaptation-20260918.md)). Gas-work sampling now
drives owner reassignment with transactional gas/particle migration and
rank-neutral continuation ([evidence](validation/spray-amr-balance-20260918.md)).
The AMR particle-interval API now handles outflow and specular reflection with
a transactional global mass/momentum/energy/species boundary ledger
([evidence](validation/spray-amr-boundaries-20260918.md)). The reacting wrapper
now combines gas and particle physical-boundary inventories transactionally
for common or recursive gas steps
([evidence](validation/spray-amr-reactive-boundaries-20260918.md)). Fixed/selected
applications now connect physical gas/particle boundaries, reservoir timestep
bounds, boundary-aware migration/regridding and cumulative boundary inventory
through changed-rank checkpoint continuation
([evidence](validation/spray-amr-open-application-20260918.md)). A stationary
wall-film migration API now preserves mixed-film payloads across finest-cell
ownership changes and regridding, including exact upper-wall positions
([evidence](validation/spray-amr-film-migration-20260918.md)). A planar mixed-film
interval API now couples heating/evaporation to owner-local AMR gas with a
transactional wall-energy ledger
([evidence](validation/spray-amr-film-step-20260918.md)). The reacting wrapper
now surrounds particle/gas stepping with film half-intervals and publishes
film state and wall heat only with whole-step success
([evidence](validation/spray-amr-film-reactive-20260918.md)). Coupled topology
remapping now carries stationary films transactionally with gas and particles
([evidence](validation/spray-amr-film-remap-20260918.md)). Gas-work-based
rebalancing also carries film state atomically with gas and particles
([evidence](validation/spray-amr-film-balance-20260918.md)). AMR particle intervals
now connect the existing planar impact model to film creation, wall impulse/energy
and optional splash children, and the reacting wrapper carries those films
([evidence](validation/spray-amr-impact-20260918.md)). Film-aware checkpoint
records now restore gas/particles/films and cumulative budgets together across
rank-count changes ([evidence](validation/spray-amr-film-checkpoint-20260918.md)).
Fixed/selected applications now expose opt-in wall impact/film settings,
film CSVs, coupled adaptation/balancing and changed-rank continuation
([evidence](validation/spray-amr-film-application-20260918.md)). A serial particle
segment query now finds first contact with the original nodal level set using
the gas EB builder's tetrahedral subdivision
([evidence](validation/spray-eb-intersection-20260918.md)). A force-free serial
population drift now reflects specular EB contacts, follows the remaining
flight time and publishes particle/wall inventories transactionally
([evidence](validation/spray-eb-drift-20260918.md)). Serial single-grid EB field
sampling and integrated gas deposition now use connected fluid-volume-weighted
stencils and transactional EOS recovery
([evidence](validation/spray-eb-field-20260918.md)). A serial EB mixture interval
now connects drag/heating/evaporation, representative-drop breakup, specular
drift and conservative gas feedback with whole-interval rollback
([evidence](validation/spray-eb-step-20260918.md)). Particle half-intervals now
wrap the existing serial reactive EB gas update, preserving separate exterior,
wall and chemistry ledgers and atomic whole-step publication
([evidence](validation/spray-eb-reactive-20260918.md)). Fixed/selected
[serial EB spray applications](../cases/eb_spray_3d/README.md) now connect initial
injection, shared nodal geometry, coupled timestep control/retries, gas/particle
CSVs and cumulative inventories
([evidence](validation/spray-eb-application-20260918.md)). Ordered timed cone
events now support tied times, zero dose, per-event mass scaling and an injected
inventory ledger ([evidence](validation/spray-eb-schedule-20260918.md)). Immutable
serial EB spray snapshots now resume gas, particles, scheduled injection and
cumulative inventories with numerical/geometry context checks
([evidence](validation/spray-eb-restart-20260918.md)). A bounded MPI ownership
primitive now migrates AMR particles with collective nodal EB fluid checks
([evidence](validation/spray-mpi-eb-tree-20260918.md)). Composite EB particle
stencils now combine exposed coarse/fine support, physical fluid volumes and
centroid visibility checks
([evidence](validation/spray-eb-tree-stencil-20260918.md)). MPI sampling now
uses these stencils to fetch only requested owner-local field values, with
collective geometry checks and transactional sample publication
([evidence](validation/spray-mpi-eb-sample-20260918.md)). MPI EB deposition now
routes integrated particle feedback by owner, restricts hidden parents with
fluid-volume weights and recovers active-cell EOS temperatures transactionally
([evidence](validation/spray-mpi-eb-deposit-20260918.md)). A distributed EB
particle interval now connects mixed-liquid source evolution, specular drift,
conservative gas feedback and final owner migration with whole-interval rollback
([evidence](validation/spray-mpi-eb-step-20260918.md)). A distributed EB gas/spray
split now combines particle half-steps with MPI reaction/transport/hydro and
commits state and global ledgers atomically
([evidence](validation/spray-mpi-eb-reactive-20260918.md)). Fixed/selected MPI EB
spray executables now connect initial mixed-liquid injection, coupled stepping,
combined budgets and per-rank particle output
([evidence](validation/spray-mpi-eb-application-20260918.md)). Scheduled injection
now supports tied, zero-scale and final-time events with an explicit injection
budget ([evidence](validation/spray-mpi-eb-schedule-20260918.md)).
A particle checkpoint payload now preserves complete parcel/component and
injection metadata across MPI rank-count changes, with rollback on rejection
([evidence](validation/spray-mpi-eb-particle-checkpoint-20260918.md)).
Combined gas/particle application checkpoints now restore on a changed MPI
rank count with pending injection events and cumulative ledgers preserved
([evidence](validation/spray-mpi-eb-restart-20260918.md)).
The application also redistributes whole patches using measured gas work and
commits gas/particle ownership together
([evidence](validation/spray-mpi-eb-balance-20260918.md)).
A coupled remap library now moves refined patches, coarsens and refines with
transactional gas/geometry/particle publication and retained original nodal
geometry ([evidence](validation/spray-mpi-eb-remap-20260918.md)).
Prescribed EB spray regrids now connect to the application and changed-rank
checkpoint continuation across stored hierarchy stages
([evidence](validation/spray-mpi-eb-regrid-application-20260918.md)).
MPI EB injection now supports independent port geometry, count, mass, diameter,
temperature and composition selected per event, with restart-bound inputs
([evidence](validation/spray-mpi-eb-ports-20260918.md)). Per-port tabulated rates
now integrate explicit packet intervals with restart-bound profile inputs
([evidence](validation/spray-mpi-eb-profile-20260918.md)). External curve-file
loading and physical qualification remain open. Liquid wall-impact support is described in the subsequent evidence. Lateral film transport remains a future extension beyond the pinned stationary-film path ([source audit](validation/film-reference-scope-20260918.md)).
Coarse/fine interpolation currently uses piecewise-constant
cell values; higher-order interface treatment/convergence remain unqualified.

Particle MPI groundwork now includes Cartesian cell ownership and stable send-buffer
partitioning with component/KH state preservation
([evidence](validation/spray-partition-20260916.md)). MPI migration and receive-side ID ordering now pass 2/4-rank checks
([evidence](validation/spray-mpi-migration-20260916.md)). Distributed gas/film
feedback and application integration remain open. Collective 64-bit ID reservation
now rejects overflow and inconsistent rank histories
([evidence](validation/spray-mpi-ids-20260916.md)). Collective computational splitting
now uses these IDs while preserving component and KH inventories
([evidence](validation/spray-mpi-split-20260916.md)). Collective cone injection now
reserves IDs and retains mixed component budgets with rollback
([evidence](validation/spray-mpi-injection-20260916.md)). Distributed field sampling
now requests only interpolation cells from owner ranks and matches the serial
stencil ([evidence](validation/spray-mpi-sampling-20260916.md)). Owner-routed
feedback deposition now matches serial stencil accumulation
([evidence](validation/spray-mpi-deposition-20260916.md)). Gas feedback now recovers
owner-block temperatures and commits conserved state and temperature collectively
([evidence](validation/spray-mpi-gas-20260916.md)). Periodic mixed-particle drag/heating/evaporation now uses owner sampling, a
common substep, collective feedback retry and end-of-interval migration with
whole-interval rollback ([evidence](validation/spray-mpi-coupling-20260916.md)).
Representative ETAB and KH-RT breakup now run inside the MPI mixed-particle
candidate update with component-inventory preservation
([evidence](validation/spray-mpi-breakup-20260917.md)). Breakup CLI settings and
deformation/RT checkpoint state now support changed-rank continuation
([evidence](validation/spray-mpi-breakup-app-20260917.md)). A collective KH shedding
API now constructs children with unique IDs, inherited composition and a gas-energy
ledger ([evidence](validation/spray-mpi-kh-release-20260917.md)). MPI coupled
feedback, input and shedding/lineage restart are now connected, including ID
rollback on rejected updates
([evidence](validation/spray-mpi-kh-coupled-20260917.md)). Periodic x-slab SGS now
shares the serial Smagorinsky face kernel, with MPI SSPRK2, coupled half-steps,
CLI coefficients and restart identity
([evidence](validation/spray-mpi-les-20260917.md)). Stationary mixed wall-film
patches can now migrate to Cartesian owners while retaining geometry and component
inventories ([evidence](validation/spray-mpi-film-migration-20260917.md)). An MPI
film-coupling API now connects heating/evaporation to owner-routed gas feedback
with a global wall-energy ledger and whole-interval rollback
([evidence](validation/spray-mpi-film-coupled-20260917.md)). MPI checkpoints now
collect and redistribute film/component state and the wall budget across rank
changes ([evidence](validation/spray-mpi-film-restart-20260917.md)). A collective
mixed-wall impact event API now reserves splash IDs and commits parent, child,
film and wall-budget updates transactionally
([evidence](validation/spray-mpi-impact-20260917.md)). Distributed film-footprint
queries now sum overlapping film heights across owners
([evidence](validation/spray-mpi-film-query-20260917.md)). An MPI planar-crossing
wrapper now computes mixture impact conditions, queries film thickness and resolves
wall events ([evidence](validation/spray-mpi-wall-crossing-20260917.md)); integration
with nonperiodic gas/particle stepping and the application remains open.
A two-layer MPI ghost API now combines physical wall/inlet/outlet conditions with
neighbor exchange, including transverse corners
([evidence](validation/spray-mpi-boundary-ghost-20260917.md)). MPI hydro now uses
these boundaries in both SSPRK2 stages and returns a global physical-face flux
ledger ([evidence](validation/spray-mpi-boundary-hydro-20260917.md)).
MPI molecular transport now shares the serial physical-boundary flux kernel,
including isothermal conduction and tangential wall work, with a global boundary
ledger ([evidence](validation/spray-mpi-boundary-transport-20260917.md)).
MPI SGS now uses the same physical-boundary ghost and wall-heat rules and
returns a global surface ledger ([evidence](validation/spray-mpi-boundary-les-20260917.md)).
The MPI particle-coupling API now supports open-face removal, stationary wall
reflection and nonperiodic owner migration with a replicated boundary ledger
([evidence](validation/spray-mpi-particle-boundary-20260917.md)). Mixed wall-impact/film
creation now participates in the coupled driver transaction
([evidence](validation/spray-mpi-impact-coupled-20260917.md)). Application impact/film
inputs, per-rank film outputs, total inventories and splash-lineage restart now
connect these models ([evidence](validation/spray-mpi-film-app-20260917.md));
physical qualification remains open.
Physical boundary data now pass through the complete MPI particle/chemistry/
transport/SGS/hydro split, with a summed boundary ledger and whole-step rollback
([evidence](validation/spray-mpi-boundary-split-20260917.md)). Application inputs now
expose physical faces and wall data; history and rank-neutral restart retain the
cumulative boundary budget ([evidence](validation/spray-mpi-boundary-app-20260917.md)).
Existing distributed films can now heat/evaporate in symmetric half-steps around
the complete coupled update, with film/gas rollback and wall-budget accumulation
([evidence](validation/spray-mpi-film-split-20260917.md)).
Periodic particle/reaction splitting now connects cell-local chemistry between
particle half-steps with collective whole-interval rollback
([evidence](validation/spray-mpi-reaction-20260916.md)). Distributed wall/breakup
events remain open. Periodic
x-slab gas advection now exchanges two halo layers and uses characteristic PLM /
SSPRK2, including optional integration between reaction half-steps
([evidence](validation/spray-mpi-hydro-20260917.md)). Molecular transport now
exchanges gas and species-limiter halos and connects as half-steps around hydro
([evidence](validation/spray-mpi-transport-20260917.md)). The fixed-mechanism
`pelef_mpi_spray_3d` executable now connects initial two-component cone injection,
global time-step selection/retry, coupled stepping, budget history and rank-local
CSV output ([evidence](validation/spray-mpi-application-20260917.md)). Scheduled
independent multiport injection now merges up to 16 port schedules with deterministic
simultaneous-event order and changed-rank restart
([evidence](validation/spray-mpi-multiport-20260917.md)). Independent per-port rate
tables now override the shared profile with complete restart identity
([evidence](validation/spray-mpi-port-profiles-20260917.md)). Finite-radius nozzle,
solid/hollow/swirl and size-distribution controls now reach the MPI injector
([evidence](validation/spray-mpi-nozzle-20260917.md)). Shared piecewise-linear injection
profiles now provide per-event integrated mass, including zero-dose event handling
and restart identity ([evidence](validation/spray-mpi-profile-20260917.md)). A shared MPI application
now also serves `pelef_mpi_spray_3d_selected`, with 21-species FFCM1_Red execution
and restart evidence ([validation](validation/spray-selected-mpi-20260917.md)).
The selected case also exposed and corrected roundoff-negative final-species mass
closure ([correction](validation/reactive-zero-species-closure-20260917.md)). Repeated timed cone events now clip time steps, retain injected budgets and
resume from the saved event count ([evidence](validation/spray-mpi-schedule-20260917.md)).
MPI spray CLI save/restart
now preserves model identity and budgets, with focused same-rank and 2↔4-rank
continuation evidence ([validation](validation/spray-mpi-restart-application-20260917.md)).
The checkpoint API now collects mixed particles and owner gas blocks into a root
file and redistributes them with a different communicator size
([evidence](validation/spray-mpi-checkpoint-20260917.md)).
The checkpoint API now preserves a 64-bit ID high-water mark with strict validation,
ready for MPI reservation state ([evidence](validation/spray-id-checkpoint-20260917.md)).

The spray CVODE branch now builds after removal of misplaced wall-film setup;
focused native/CVODE chemistry tests pass
([repair evidence](validation/spray-cvode-build-repair-20260916.md)).

## Repository consolidation implemented

The README is an entry point rather than a chronological implementation log.
Contributor instructions, AGENTS.md, the documentation index and PR template
separate current capability, historical evidence and future work.

The 2613-line root CMake and 9138-line test registration were split into five
root target modules and nineteen test modules, preserving include order and
directory scope. Generated target graphs, install rules and CTest commands
matched in seven clean configurations; see [the layout record](validation/build-layout-20260909.md)
and [CMake navigation](../cmake/README.md).

Hosted workflows cover source provenance/regeneration, complete serial and MPI
suites, selected 0D CVODE with pinned SUNDIALS 7.2.0 and Cantera 3.2.0, and a
separate tests-disabled Release installation. Source identity, dependency
versions, JUnit and logs are retained. An implemented workflow is not itself a
passed qualification; completed results must be associated with its tested SHA.

Historical draft #11 was reviewed and closed as superseded: its 17 numerical,
mechanism, application/case and verification paths were byte-identical at the
later integrated full-chemistry commit. Its unique decision note is preserved
[in the historical record](history/pr11-implicit-chemistry-decision.md).
No historical branch was deleted. The audited ancestor-only bulk cleanup was
blocked by the connection's safety restriction, not performed silently.

## Evidence ledger

The saved checkpoint is `509d07929d52d9e570a6af1f6ec7e5836bd4fe5d`.
Initial PR CI failed because its vendored YAML description lacked two spaces.
Restoring the exact Cantera package bytes repaired source provenance without
changing the mechanism bundle or generated Fortran; see [the integrity repair](validation/checkpoint-integrity-20260908.md).

A later targeted numerical review found a false-success return in NASA7 molar
conversion overflow guards. Commit `bb2980e52c2e6a4bfed353dd280fe1f0e21f4893`
corrects the failure flag without changing physical formulas or tolerances.
New cp/enthalpy/entropy rejection cases fail on the parent implementation and
pass after the fix; seven focused GNU14 Debug tests passed. See [red/green evidence](validation/nasa7-conversion-status-20260909.md).
This correction is intentionally separate from the behavior-preserving layout.

**Reported historical local qualification:** the [2026-09-08 handoff](handoff-pelef-0.245.0-20260908.md)
reports 4570/4570 across Debug 524, Release 524, Debug+Cantera 529, MPI Debug 694,
MPI Release 694, CVODE Debug 533, CVODE Debug+Cantera 539, and CVODE Release 533.
It also reports tests-disabled 42/42 ELF, separate tests-enabled 23/23 installed
restart checks, and installed application/checker/byte comparisons of 32/32,
8/8 and 12/12. These belong to that frozen worktree and its local evidence;
the raw server directories were not retrieved here. Interrupted 227-test and
excluded literal-prefix exploratory results are not successful qualification.

**Independent source-integrity qualification:** commit
`566010ffc77a6de25d452d084ff61fa2bb533d95` passed
[run 34208154353](https://github.com/GTT10/Fortranslate/actions/runs/34208154353),
including 14 checker regressions, exact hashing and full bundle regeneration.
Final serial/MPI/CVODE/clean-install results belong to PR #102 and their own
source identities. They must not be replaced by older source successes.

The earlier focused 10/10 serial and 23/23 MPI results remain in
[validation/0.245.0.md](validation/0.245.0.md). Historical evidence is preserved,
not overwritten or added together as though all configurations were identical.

See [regular 3D boundary inputs](../cases/reactive_boundaries_3d/README.md) for the new
hydro/chemistry boundary scope and flux accounting.

## Remaining work

[The roadmap](completion_roadmap.md) distinguishes further solver development
from repository consolidation. General 3D AMR/EB, production detailed-fuel stiff CFD,
full LES/spray qualification and external physical validation remain genuine engineering work,
not features implicitly completed by housekeeping or large test counts.
