# Experimental regular 3D spray combustion (0.246.0)

This is a working opt-in Eulerian--Lagrangian path, not production diesel-spray
qualification. It couples the existing regular 3D gas solver to
single-component point parcels, molecular transport, a deviatoric Smagorinsky
closure and native or CVODE cell chemistry. It does not replace the existing
AMR/EB/MPI applications or claim those features are connected to particles.

## Particle source export for SDC integration

MPI AMR/EB deposition and `advance_mpi_amr_mixture_coupling` accept optional
`gas_change`: an owner-local field of accepted, time-integrated gas feedback per
fluid volume. It is accumulated independently of the gas background and published
with the accepted gas/particle update. Failures retain its previous value. Divide
by the interval for a source rate; wall/outflow inventories remain separate.
The EB SDC gas driver now accepts old/new source callbacks and reports their
inventory separately from external forcing. The common-step MPI EB particle
wrapper now connects kick/drift stages, physical gas boundaries and final
computational splitting. The MPI EB spray application exposes common-step and recursive SDC. The new-time particle half kick runs only at the final SDC iteration; intermediate iterations retain the old particle source ([reference-order correction](validation/spray-sdc-reference-order-20260927.md)). See the
[source-export checks](validation/spray-sdc-source-export-20260927.md) and
[gas-stage checks](validation/sdc-eb-auxiliary-source-20260927.md).

The library additionally provides `advance_mpi_eb_spray_sdc_subcycle`, connecting
recursive gas stages, physical/virtual/ghost particles, predictor halos and final
physical collection for a nonperiodic drag/heating/evaporation core. The MPI spray
application selects this path with `time_integrator='sdc'` and `subcycling=.true.`
([application checks](validation/spray-recursive-app-20260927.md)). Recursive ETAB and KHRT diameter evolution now uses the existing mixture
breakup kernels on all particle populations. Computational splitting runs after
outer synchronization, with collective ID allocation and rollback
([breakup checks](validation/spray-recursive-breakup-20260927.md)). KHRT
shedding children now run on old-time physical stages, with their transverse
energy charged to the gas and IDs committed with the accepted recursive step
([child checks](validation/spray-recursive-children-20260927.md)). The interior
child placement formula now has a direct pinned-source and 2/4-rank stage check
([placement audit](validation/spray-child-placement-20260927.md)); child count
and angular-distribution parity remain open. Recursive films remain open; their
frozen stage-time halo exchange kernel now returns same-level and physical
source inventories transactionally
([kernel checks](validation/spray-film-halo-20260927.md)). Uniform particle acceleration now
brackets each level interval with old/new half kicks; only physical particles
contribute to the external force ledger
([acceleration checks](validation/spray-recursive-acceleration-20260927.md)). End-of-level virtual cleanup, final-subcycle ghost cleanup and
synchronized physical subtree redistribution now run after gas synchronization
([lifecycle checks](validation/spray-recursive-lifecycle-20260927.md)). Inflow metadata and isothermal
wall temperatures now reach recursive gas stages; particle faces map inflow to
outflow and isothermal walls to no-slip reflection, as in the common-step wrapper.
See the [coupled core checks](validation/spray-recursive-core-20260927.md) and
[boundary checks](validation/spray-recursive-boundary-20260927.md).

The recursive population source uses the pinned PelePhysics containing-cell
rule; gas sampling still uses surrounding donors. At outer synchronization,
PeleF replaces the provisional copy-supported leaf exchange with the actual
physical-particle exchange, deposited at its containing composite leaf cell.
This conservative synchronization is an addition to the pinned same-level
transfer, not a claim of exact upstream numerical parity. Both directions of
coarse/fine crossing now pass the unchanged inventory criterion; see the
[source synchronization checks](validation/spray-recursive-source-sync-20260927.md).
The preceding [source audit and failure record](validation/spray-recursive-crossing-20260927.md)
remains historical evidence. Three-level coupled evolution, temporal accuracy,
strong-coupling admissibility and the remaining recursive particle models remain
open. A two-level application case now verifies scheduled regridding followed
by 2-to-4-rank restart, reaction-source history and a pending injection event
([continuation checks](validation/spray-recursive-restart-regrid-20260927.md)).
The MPI spray application's load-balancing path also transfers SDC reaction
history with gas and particle ownership, including rollback after particle or
film migration failure ([focused checks](validation/spray-sdc-rebalance-history-20260927.md)).

## Level-bounded particle ownership

`locate_amr_tree_spray_cells`, `partition_amr_tree_spray_parcels`, and the
regular/EB MPI migration APIs accept optional `maximum_level` (zero based).
It selects the finest covering node no deeper than that level, including
parents covered by refinement. Omission retains finest-level ownership.
MPI ranks must agree on the effective limit; invalid limits or disagreement
leave parcels, components and partition outputs unchanged.

This permits a caller-owned particle copy to be assigned to coarse-level
owners without coarsening the gas topology or changing parcel coordinates,
IDs, mass or breakup history. Virtual population creation is described below;
ghost creation is described below. Aggregation and the connection to recursive
SDC remain required. The pinned reference's
`Source/Particle.cpp::particleMKD` creates virtual/ghost representations before
level updates; these ownership APIs provide a prerequisite for that lifecycle.
See the [focused checks](validation/spray-level-ownership-20260927.md).

## Coarse virtual particle populations

`create_mpi_amr_level_particles` extracts the physical population whose finest
containing node is exactly `target_level`. Unlike coarse virtual copies, these
level populations partition the input IDs across levels. It returns owner-local
parcels, components and matching node indices, with all outputs replaced only
on success. The original population stays read-only; these snapshots must not
be added back as extra physical particles. Optional EB validation uses the
existing migration path. This supplies level-specific physical inputs for the
recursive driver; their clocks and final merging are still to be connected.
See [physical-level checks](validation/spray-physical-level-20260927.md).

`spray_level_population` groups a level's particle/component arrays and node
indices. `collect_mpi_amr_level_particles` combines synchronized physical
level snapshots, removes dead entries and migrates live particles according
to their updated positions. Stale node indices are ignored. Before geometric
migration, an ID-based exchange detects duplicate live IDs globally, even at
different positions. Any failure preserves the caller's full population and
components; inputs remain read-only. Each local input and temporary ID-owner
bin has the existing one-million-particle bound. Ghost/virtual copies must
not be supplied. Synchronization timing and invoking collection from the
recursive driver remain caller responsibilities.
See [collection checks](validation/spray-level-collection-20260927.md).

`create_mpi_amr_virtual_particles` constructs a separate, nonaggregated
population for a requested `target_level`. It selects physical parcels whose
finest covering node is deeper than the target and migrates copies to the
coarse owners. Coordinates, component inventories and physical/breakup state
are copied, while the source population remains unchanged. Supplying
`nodal_phi` uses the EB migration and solid-point validation path.

Outputs replace the previous virtual population only after collective success.
The finest level produces an empty population; target disagreement, invalid
inputs, capacity exhaustion and EB rejection retain existing outputs.
The input must contain each active physical parcel once, with sorted local IDs.
Source IDs are retained for correspondence; virtuals must stay separate from
physical inventories, injection and checkpoints. Callers must regenerate/clear
these temporary populations at the level lifecycle boundaries.

This implements the unaggregated coarse representation needed by the pinned
PeleC virtual-particle lifecycle. It does not implement AMReX cell aggregation,
reserved virtual IDs, level-specific interpolation/deposition or recursive time
integration. Ghost creation is described below. Those connections remain open.
See [creation checks](validation/spray-virtual-particles-20260927.md).

## Fine-level ghost particle populations

`create_mpi_amr_ghost_particles` copies immediate-coarse physical parcels
into fine-level boxes enlarged by `grow` fine cells. `target_level` must be
at least one. The separate output contains parcel copies, component masses
and target node indices. Periodic images participate in lookup without
changing parcel coordinates; valid boxes take precedence over enlarged boxes,
and overlapping enlarged boxes use tree order for a single destination.
`locate_amr_tree_ghost_cells` exposes the same lookup, with node zero for no
match and local cell indices that may lie outside the valid patch.

The source is read-only. Creation first validates a private physical migration,
then selects only parcels whose finest covering level is `target_level-1`.
Optional `nodal_phi` retains EB solid-point validation. All output arrays
publish together after successful exchange; capacity failure, invalid geometry
or rank-disagreeing target/width retains previous outputs. Source IDs are
retained, so ghosts must remain separate from physical inventories and
checkpoints. These copies require level-local ghost gas fields and source
transfer; they cannot be fed directly into the existing composite coupling
and treated as additional physical particles.

This supplies the creation step used by the pinned reference's
`setupGhostParticles`; grown-box tie ordering is this implementation's explicit
contract and has not been compared to AMReX locator ordering. Per-subcycle
refresh/removal, ghost evolution, source transfer and recursive application
routing remain open. See [creation checks](validation/spray-ghost-particles-20260927.md).

`redistribute_mpi_amr_ghost_particles` reassigns an existing ghost population
after movement or an owner-map change. It retains live copies in the target
level's valid or grown boxes, including copies that have entered fine valid
cells. Dead entries and copies outside every grown box are removed. This is
temporary-copy retirement, not loss of physical particle inventory. Returned
node indices follow the exchanged, ID-sorted population. Capacity or validation
failure retains particles, components and node outputs together.

The operation reuses the creation path's private migration checks; optional
nodal EB geometry still excludes periodic EB domains and rejects solid live
particles. It does not create additional copies from the physical population
or choose the recursive lifecycle at which regeneration should occur.
See [redistribution checks](validation/spray-ghost-redistribution-20260927.md).

`advance_mpi_amr_ghost_halo` connects one ghost-population stage: halo-based
kick, drift, same-level reaction summation into valid cells, then ghost
redistribution. Kick and drift durations are explicit. The caller supplies
stage-time gas halos, target level, grown-box width and ghost capacity. Updated
particles, components, node indices, source/remainder fields and boundary
ledgers publish together only after final redistribution succeeds.

The ghost stage disables direct coarse routing. Contributions without a
same-level destination remain in `remainder`, available for level-local halo
source use. They must not be added as a second coarse physical-particle source.
This follows the pinned PeleC/PelePhysics transfer path, which sums same-level
boundaries and copies into the current level's source, including requested
ghost cells. The generic halo kick's explicit coarse-routing option remains
available for callers with separately established accounting. This corrects
the earlier integrated-stage coarse-routing behavior; see
[source-scope correction](validation/spray-ghost-source-scope-20260927.md).

This stages an existing ghost population; recursive clocks, creation/removal
at lifecycle boundaries, breakup/children and integration of ghost feedback
without duplicate physical/virtual contributions remain to be connected.
See [stage checks](validation/spray-ghost-stage-20260927.md).

`advance_mpi_eb_ghost_stage` additionally fills the stage-time conserved and
temperature halos directly from the gas state. It forwards optional ancestor
state/temperature endpoints and interpolation fractions to the existing gas
halo implementation, with snapshot restriction disabled. Caller-supplied
target-resolution volume fractions and centroids provide EB sampling metrics.
It then runs the complete ghost stage above. Gas snapshots remain read-only;
halo-fill failure leaves all particle/source/ledger outputs untouched.
The underlying EB redistribution still requires nonperiodic geometry.
See [time-interpolated stage checks](validation/spray-ghost-temporal-20260927.md).

The MPI recursive gas driver accepts optional `halo_source` using the
`mpi_eb_sdc_recursive_halo_source` interface. It is an alternative to
`auxiliary_source` and uses the same required `auxiliary_change` ledger and
optional cross-level routing. In addition to interval start time, stage and
depth, it supplies read-only ancestor state/temperature endpoints and the
fractions needed by the gas halo fill. Stage zero receives interval-start
fractions; positive stages receive interval-end fractions. Fractions for
nonancestor levels are zero. MPI ranks must agree on the callback kind.
The callback still owns private particle work until the outer gas step succeeds;
this interface alone does not connect the recursive particle lifecycle.
See [callback checks](validation/sdc-particle-halo-context-20260927.md).

`compute_sdc_mpi_eb_tree_hydro_3d` accepts optional `predictor_halo` for a
complete, caller-assembled three-cell predictor-source halo. MPI ranks must
agree on its presence. Owner allocation, shape, finiteness and exact agreement
with the normal predictor on valid fluid cells are checked before use. Only
ghost values may differ; the supplied field replaces the normal halo rather
than being added to it. Packed or physical lower bounds are accepted.
This lets a particle-aware caller retain directly deposited halo sources in
CTU prediction. Assembling that complete predictor and passing it from the
recursive particle callback still require integration.
See [predictor halo checks](validation/sdc-predictor-halo-20260927.md).

## Level-bounded gas sampling

Regular/EB cloud stencils, the serial regular sampler and MPI regular/EB
sampler accept optional `maximum_level`. The cloud width uses the containing
cell at that level ceiling; deeper nodes are excluded and their covered
parents become available. MPI exchanges only the resulting requested cells
and requires all ranks to agree on the effective ceiling. Omission retains
the composite finest-grid interpolation. Invalid ceilings and failed queries
retain the previous sample output, including empty-query cases.

This allows coarse virtual populations to read the corresponding gas view.
It is a composite view capped at a level, not a fine-level halo interpolator:
coarse regions outside that level's patches still use coarse cells. Time
interpolation of parent endpoints, fine ghost gas reconstruction, and separate
level source deposition/transfer remain required for recursive spray SDC.
The gas/EOS deposition API and particle evolution now expose the same ceiling
for provisional level views, as described below.
See [sampling checks](validation/spray-level-sampling-20260927.md).

## Level-bounded particle source deposition

`deposit_mpi_amr_tree_spray_field` accepts optional `maximum_level` for
temporary source fields. An explicit ceiling uses the same capped cloud as
sampling and adds volume-normalized contributions only to its exposed cells.
It skips parent restriction so deeper fields cannot overwrite source increments
stored in covered coarse cells. `gas_change` exports these raw increments
independently of the background field. Fine fields above the ceiling remain
unchanged. MPI agrees on both ceiling presence and value before deposition.

The optional ceiling selects source-field semantics, even when it equals the
finest level: covered ancestors are not refreshed. Omitting it retains physical
composite deposition and parent restriction. This supplies per-level source construction;
source transfer/synchronization and the recursive spray driver still need to
consume the separate level fields without double-counting virtual/ghost copies.
See [source checks](validation/spray-level-source-20260927.md).

## Evolution in a capped gas view

`advance_mpi_amr_mixture_coupling` accepts optional `maximum_level`. Sampling,
gas/EOS deposition, the spatial substep limit and final parcel migration use
the same ceiling. `deposit_mpi_amr_tree_spray_gas` also exposes this option;
it checks the gas EOS transactionally and preserves finer temperatures exactly.
The returned gas is a provisional capped view: no parent restriction is done,
even for an explicit finest-level ceiling. `gas_change` exports its accepted
increments. Omission retains the original composite coupling.

This permits evolution of separate coarse virtual populations. The caller
must provide the intended population and synchronize sources at the appropriate
time; the routine does not select or deduplicate physical/virtual particles.
It still uses capped composite interpolation, so fine ghost populations require
the separate halo operators to be connected. Wall-film populations now use
the same ownership ceiling, as described below.
See [evolution checks](validation/spray-level-evolution-20260927.md).

## Liquid films in a capped gas view

`migrate_mpi_amr_wall_films` accepts optional `maximum_level`, selecting the
finest containing cell at or below that ceiling. Exact wall coordinates,
normals, material and component data are preserved. EB wall intersection
validation still uses the original nodal geometry. Ranks must agree on the
effective ceiling; invalid ceilings or failed exchange retain both outputs.

`advance_mpi_amr_film_coupling` forwards the ceiling to migration, gas sampling
and deposition, and uses the capped spacing for planar wall gas probes.
`advance_mpi_amr_mixture_coupling` also forwards it when migrating impact films.
The gas result is a provisional capped view awaiting synchronization, with
separate `gas_change` and wall-energy inventories. Recursive film scheduling
and physical/virtual population accounting remain caller responsibilities.
See [film checks](validation/spray-level-film-20260927.md).

## Particle sampling from complete gas halos

`kick_mpi_amr_mixture_halo` connects halo sampling to the existing mixed-liquid
drag, heating and evaporation kernel for one caller-selected kick interval.
It accepts assigned node indices, conserved and temperature halos, and optional
fluid-fraction/EB visibility data. Gas remains read-only; the output is a fresh
source halo per fluid volume using the same weights as sampling. Particle,
component and source outputs publish together after collective success.

Without optional drift, positions remain fixed. IDs and assigned nodes are
retained in both modes. The caller supplies stage-time halos, chooses the kick
duration (including any half-step factor), and handles breakup/children, dead-particle removal and source
transfer. This is a kick building block for physical, virtual or ghost
populations, not a complete recursive step. Source application/EOS validation
and retries belong to the enclosing gas/particle transaction.
See [halo kick checks](validation/spray-halo-kick-20260927.md).

Supplying `drift_interval`, `boundary_faces`, `wall_change` and `boundary_change`
together adds a subsequent drift with the kicked velocity. The kick and drift
durations are independent (for example, half an interval and one interval).
The existing level-set drift handles EB specular contacts and physical
periodic, reflecting and open boundaries; regular geometry uses an all-fluid
level set. EB and outer-boundary inventories are returned separately, including
the vapor-species inventory removed by outflow. Dead entries remain with zero
component mass. Reaction deposition uses the original sampling positions.
Drift failure rolls back the preceding kick, all source fields and both ledgers.
MPI ranks must agree on the enabled arguments, duration and boundary labels.
This does not migrate particles to new owners or regenerate ghost populations;
the caller must keep assigned nodes and halo support valid for subsequent kicks.
See [drift checks](validation/spray-halo-drift-20260927.md).

Supplying `target_level`, `level_source` and `source_remainder` together also
transfers the raw reaction halo before publication. Every queried node must
belong to that level. The returned `level_source` is a fresh owner-local field
on valid cells; `source_remainder` contains untransferred halo contributions.
Optional `route_to_coarse` enables conservative transfer to covering ancestors,
and `maximum_received` bounds the exchange. EB transfer requires `eb_geometry`
together with the sampling `fluid_fraction`, preserving fluid-volume scaling.
The raw `source_halo` is still returned for diagnostics. Transfer failure rolls
back particles, components, all source outputs and optional drift ledgers.
Applying the source to gas, synchronizing levels and avoiding duplicate
physical/virtual/ghost contributions remain driver responsibilities.
See [integrated transfer checks](validation/spray-halo-transfer-20260927.md).

`sample_mpi_amr_particle_halo` reads a caller-filled owner-local halo at the
particle's assigned node resolution. Each query supplies its target node;
that node must belong to the querying rank. Trilinear weights use the fine
cell centers, including halo cells outside the valid patch. Periodic queries
use the image nearest the patch center without modifying particle coordinates.
Both physical spatial bounds (`1-width:N+width`) and packed bounds are accepted.
An optional matching `fluid_fraction` halo multiplies and renormalizes weights;
zero-fluid cells are ignored, and zero total support rejects the query.

The caller must fill all required corner donors at the stage time. The existing
`fill_mpi_eb_tree_halo_3d` produces compatible layouts and accepts ancestor time
endpoints; requesting its volume-fraction output fills corners as well as faces.
Wider EB halos require its target-resolution volume metrics. This sampler does
not fill those arrays or perform time interpolation itself. Optional nodal
geometry and centroid halos enable the visibility path described below. Width/geometry disagreement,
nonfinite requested values, wrong ownership or insufficient halo support reject
collectively and preserve the previous samples. See
[halo sampling checks](validation/spray-halo-sampling-20260927.md).

## Stage-time EB gas sampling for particles

`sample_mpi_eb_tree_particle_gas` connects the existing MPI EB gas halo fill to
particle-site sampling. It accepts ancestor state/temperature endpoints and
fractions, fills full halos including corners, samples conserved variables and
temperature guesses at assigned particle nodes, then recovers particle-site
temperature from the interpolated conserved state through the EOS. Conserved
samples and temperatures publish together only after collective success.
Parent restriction is disabled so the supplied stage states are retained.

The wrapper forwards existing boundary/inlet controls and optional target
volume-fraction halos. Wider halos retain the gas fill's requirement for
explicit target-resolution metrics. This connects stage-time parent donors to
ghost-particle gas references. Optional centroid visibility is described below;
particle evolution, source transfer and recursive application stepping remain open. See
[temporal halo checks](validation/spray-temporal-halo-20260927.md).

## EB visibility for particle halo sampling

The halo sampler accepts paired `nodal_phi` and `centroids` arguments alongside
its fluid-fraction halo. Centroids are absolute physical coordinates in an
owner-local three-component halo with the same spatial shape as the gas halo.
The nodal field must agree across ranks, match the tree's finest lattice and
satisfy periodic endpoint equality. Only donors with nonzero fluid support
and an unobstructed segment from the particle to the centroid contribute;
weights are renormalized after blocked donors are removed. Solid queries or
loss of all support reject collectively. Periodic segments reuse the existing
composite stencil's boundary-splitting intersection routine.

The stage-time wrapper forwards these arguments and requires explicit
`target_volume_fraction` when visibility is enabled. The caller must supply
centroids and fractions from the target-resolution extended geometry; the gas
halo fill retains exact interior-volume consistency checks. The sampler checks
shape and requested centroid finiteness, but cannot prove that caller-supplied
centroids came from that geometry. Without the optional pair it retains the
fraction-only interpolation path. See
[visibility checks](validation/spray-halo-visibility-20260927.md).

## Same-level particle source halo summation

`sum_mpi_amr_particle_source_halo` adds valid/halo source contributions to
matching valid cells at `target_level`, including periodic images and remote
patch owners. It uses integer cell coordinates and adds in source-rank order.
The input halo is read-only; the updated source field and `remainder` publish
together. Matched entries are zeroed in the remainder. Contributions without
a same-level receiver, including coarse/fine gaps and nonperiodic domain
exteriors, remain at their original halo location for subsequent handling.
Other levels in the source field are unchanged; no parent restriction runs.

Inputs are additive source contributions, not already synchronized ghost
copies. Every representation of the same physical cell must use the same
normalization and EB volume metrics. This routine performs no volume scaling,
wall transformation or virtual/ghost de-duplication. Optional coarse routing
is described below. The
caller must consume each halo once and handle the remainder before claiming a
completed source transfer. This supplies the same-level summation used by the
pinned PelePhysics `transferSource`/`SumBoundary` path; recursive integration
and synchronization scheduling remain open. See
[source halo checks](validation/spray-source-halo-20260927.md).

## Conservative coarse routing of particle sources

`sum_mpi_amr_particle_source_halo` additionally accepts `route_to_coarse`
(default false). If no same-level cell receives an in-domain contribution,
true selects the finest covering ancestor and scales the source density by
source-cell volume / destination-cell volume. Optional paired `eb_geometry`
and `halo_volume_fraction` replace geometric volumes by fluid volumes. The
replicated EB hierarchy is validated; source halo fractions are owner-local
and must match their target-resolution geometry. A routed nonzero source with
zero source or destination fluid volume rejects collectively. Overflows and
zero scaling also reject before publication.

Same-level and coarse routing share the atomic source/remainder update.
Nonperiodic exterior contributions remain unresolved. This transfers source
inventories; it does not apply them to gas thermodynamics or select which
real/virtual/ghost populations contribute to each subcycle. The recursive
caller still must coordinate population clocks, source accounting and gas
synchronization. See [coarse transfer checks](validation/spray-coarse-source-20260927.md).

## Matched particle halo sampling and deposition

`deposit_mpi_amr_particle_halo` uses the halo sampler's exact normalized
weights to deposit extensive particle contributions into a fresh source halo.
Each contribution is divided by the receiving cell's fluid volume. Assigned
nodes, periodic images, fluid fractions and optional EB centroid visibility
follow the same rules as sampling; blocked donors receive no source. It uses
the donor gas halo as read-only input and publishes the source halo only after
collective success. Previously populated outputs survive a rejected trial.

The sampling API also accepts paired `contributions` and `source_halo` options
for joint sampling/deposition; both outputs publish together. The deposition
wrapper calls this shared operation. Resulting source halos can be passed to
`sum_mpi_amr_particle_source_halo`, using the same volume metrics, for MPI and
coarse routing. This does not evolve a parcel or choose whether a physical,
virtual or ghost contribution is included in a level's inventory. Those
population and subcycle decisions remain the recursive driver's responsibility.
See [matched-weight checks](validation/spray-halo-deposition-20260927.md).

## AMR and embedded boundaries

An AMR ownership library now locates parcels in the finest covering regular
3D patch-tree cell and constructs stable rank-grouped send buffers using a
caller-supplied owner map. It handles periodic lookup and reassignment after
regridding while preserving parcel/component state. An MPI migration API now
uses this ownership and the shared parcel exchange to transfer particles to
their gas patch owners. A serial composite-field sampler/depositor is now
available, and owner-local MPI sampling now queries only the required composite
cells. MPI deposition now routes conservative contributions to those owners;
MPI gas temperature recovery now commits with the feedback. A periodic MPI AMR
mixed-particle interval now connects sampling, particle evolution, feedback and
final migration. A reacting-step wrapper also connects the existing AMR gas
solver. Fixed/selected [MPI AMR spray executables](../cases/mpi_amr_spray_3d/README.md)
now connect initial or timed injection, coupled steps, prescribed or temperature-adapted regridding and gas/particle
CSV output. Checkpoint/restart now preserves the gas/particle state and event
position across MPI rank-count changes.

The tree runtime also exports `spray_eb_intersection_mod`'s
`intersect_spray_level_set(phi,lower,upper,start,finish,hit,fraction,position,normal,cell,ok)`.
It finds the earliest fluid-to-solid segment contact using the same six
tetrahedra per Cartesian cell as `build_level_set_eb_geometry_3d`. Supply the
original nodal level set (`phi>0` is fluid), not aggregated EB face moments.
The returned unit normal points into fluid; cell indices are one-based.
Both endpoints must lie inside the closed domain and the start must be fluid
or on the surface. Surface-inward motion hits at fraction zero. A valid miss
returns fraction one, the endpoint, and zero normal/cell. Invalid input preserves
all result arguments. Geometry is bounded to one million cells, and the query
scans the segment's cell bounding box.

`spray_eb_drift_mod` adds
`drift_spray_level_set(phi,lower,upper,dt,maximum_hits,parcels,change,hits,ok)`.
It advances a population at constant velocity between stationary specular EB
contacts, reflects the normal velocity and continues the remaining time, with
a per-parcel contact limit (0 through 1024). `change(1:5)` is the represented
particle mass/momentum/energy difference; the wall receives its negative.
Liquid mass, temperature, multiplicity and component identity are unchanged.
Dead zero-mass parcels are skipped. The population, change and total hit count
publish together; failure preserves all three. Each proposed segment, including
its unreflected endpoint, must stay in the closed root domain. This API does not
combine outer-domain crossings with EB events. It remains a force-free serial
drift API; liquid impact/film models are not connected to this drift routine.

`spray_eb_impact_mod` exports `impact_mixture_level_set` for one earliest
mixed-liquid EB contact. Supply the previous position, a parent at its trial
endpoint, component masses, molecular weights, local pressure, wall temperature,
film thickness and contact angle. It uses the facet normal and the existing
planar model's Weber/Reynolds conversion and impact correlations. Successful
contacts publish the retained parent, splash children, deposited film and their
component inventories together. `wall(1:5)` is received by the wall (opposite
to the specular drift routine's system-change convention). A valid miss sets
`fraction=1`, `regime=-1`, `wall=0` and preserves population/film buffers;
invalid input preserves all state and diagnostics except `hit`/`ok`, which are
false. The caller owns film-thickness sampling, global child-ID reservation,
population insertion and remaining flight. This single-contact API is not yet
connected to EB interval/application stepping or film evolution. See
[focused impact evidence](validation/spray-eb-impact-20260918.md).

`advance_mixture_level_set_impacts` advances an allocated mixed-parcel population
through a force-free interval using this contact model. Inputs include common
pressure, wall temperature, liquid viscosity/surface tension, the film population,
an ID watermark and capacities for all intermediate parcels and films. Input
parcel IDs must be strictly increasing and no greater than the watermark.
Children receive consecutive IDs and advance for their parent's remaining time;
consumed parents are removed at publication. The contact limit applies along
each parent/child trajectory. Existing and newly deposited stationary films
contribute thickness using the planar model's normal/proximity criterion.
Processing follows parcel ID order with appended children; this is not a globally
time-ordered event queue. Capacities include consumed intermediate parents.
All population/component arrays, films, watermark, wall ledger and hit count
publish together or remain unchanged. `wall` again means received by the wall.
Pressure/properties are frozen; there is no gas feedback, film evolution, MPI
reservation or application connection in this API. Proposed flight endpoints
must remain inside the closed root domain. See
[interval evidence](validation/spray-eb-impact-interval-20260918.md).

The interval API also accepts `outflow_change(5+nc)`. Its presence enables
exterior outflow on every root face: each trajectory is clipped at its first
root crossing before EB impact dispatch. A prior EB impact retains the correct
remaining time for both parent and splash children. Escaped records are removed
alongside consumed parents, without recycling IDs. The output is negative system
change in mass, momentum, total energy and component masses (component order,
not gas-species indices); the existing `wall` output remains wall receipt.
The depth limit includes exterior events, and failures preserve all populations,
films, IDs and both ledgers. Serial source and reactive coupling now connect this
mode, including the serial application; MPI wiring remains pending. See
[impact/outflow evidence](validation/spray-eb-impact-outflow-20260918.md).

`advance_eb_mixture_coupling_allocated` now exposes the serial EB source interval
with allocated parcel/component arrays. Supply `wall_options` (the
`spray_eb_wall_options` type), allocated `films` and `film_mass`, and `last_id`
together to enable impacts. The original fixed-population
`advance_eb_mixture_coupling` remains the specular interface and calls the same
source core. Each substep samples the gas, advances drag/heating/evaporation and
representative breakup, then advances EB impacts at the updated velocity.
Gas exchange is deposited at pre-step positions. Splash children enter gas
sampling on the next substep. The optional `pressures` vector of the impact
interval carries each sampled local pressure; descendants inherit that pressure
for their remaining flight. Wall temperature and impact properties are common.
The coupling ledger is system change, hence the negative of the impact kernel's
wall receipt. Gas, particles, films, component inventories and ID watermark
commit together after the whole interval succeeds. MPI and application wiring
remain separate work. See
[coupling evidence](validation/spray-eb-impact-coupling-20260918.md).

The allocated serial source also supports KH shed-mass children with a shared
`last_id`, with or without impact films. It checks each source-updated trajectory
against the original level set before shedding: a crossing skips KH for that
substep and takes the wall path. Unhit parents may release four children through
`advance_mixture_khrt_shedding`; their transverse kinetic-energy change feeds the
parent gas stencil, and their pressure inherits the parent's sampled value.
Parents and children then drift for the source interval. IDs and pending KH
state remain private through retries and whole-interval failure. This remains
source-first operator splitting, not globally time-ordered breakup. The serial
reactive wrapper and application now pass the private shared watermark through
both source halves and retain KH state through restart; restore validates
four-child ID groups after subtracting injection and computational split IDs.
See [application evidence](validation/spray-eb-khrt-application-20260918.md).
See [serial KH source evidence](validation/spray-eb-khrt-source-20260918.md).

Set `spray_eb_wall_options%evolve` to enable stationary mixed-film heat and
evaporation in that allocated source interval. `conductivity` is a common
constant liquid conductivity (default 0.6 W/(m K)); `minimum_height` defaults to
1e-7 m and follows the existing film model's depletion rule. The independent
`spray_eb_film_step_mod` API `advance_eb_film_coupling` samples/deposits through
the EB fluid stencil, adapts its own film substeps and removes depleted films.
Film positions and normals must match the original nodal surface: a short
normal segment must cross at the patch position with the same fluid-facing
normal. Sampling is offset into fluid by the larger of 1e-8 times the minimum
cell width and 128 machine epsilons times the maximum absolute domain coordinate.
Both probe endpoints must stay inside the root domain.
The whole film interval publishes gas, film components and the wall heat ledger
together. In the combined source path, particle gas feedback precedes existing
film exchange, followed by impacts; new deposits begin evolving on the next
particle substep. Film temperature/conductivity are common settings, with no
lateral film transport, surface mesh or variable-conductivity fit in this EB API.
See [film coupling evidence](validation/spray-eb-film-coupling-20260918.md).

`advance_eb_spray_reactive_allocated` composes these allocated particle/film
half-intervals around the existing reactive EB gas step. It accepts the same
gas controls as `advance_eb_spray_reactive`, plus optional `wall_options`,
`films`, `film_mass` and `last_id` supplied together. Both source halves and
the gas update operate on private state; only a successful whole step publishes
gas, particles, films, IDs, limiter and exterior/wall/chemistry ledgers.
The existing fixed-population API delegates to this core in specular mode.
New first-half deposits participate in film exchange in the second half.
MPI support for these impact/film arguments is still open.
See [reactive impact evidence](validation/spray-eb-impact-reactive-20260918.md).

The fixed/selected serial EB spray applications now expose `wall_impacts` and
the `walls` options in `&eb_spray`, output stationary film CSVs and retain films,
component inventories, retired-parent ancestry and ID watermarks in coupled
snapshots. See the [input guide](../cases/eb_spray_3d/README.md) and
[application/restart evidence](validation/spray-eb-impact-application-20260918.md).
The checkpoint context is now V2; older V1 contexts are rejected. Existing wet
correlations can reject normal-incidence thick-film states; this boundary is
documented in the evidence and is not removed by the application connection.
The contact, impact-interval, allocated source and allocated reactive APIs expose
an optional `model_rejected` output for inadmissible classification, child-mass
partition or outgoing normal velocity. They retain transactional state rollback.
The source and application propagate this condition without timestep retry;
other step failures retain their existing retry behavior. The flag is reset
on entry and is separate from `ok`. See
[rejection evidence](validation/spray-eb-impact-rejection-20260918.md).

`migrate_mpi_amr_wall_films` now accepts optional `nodal_phi` for stationary
EB films. It validates each film against the nodal surface and assigns it to
the finest cell containing a short fluid-side normal probe; the transmitted
film retains its exact surface position. This resolves ownership at EB surfaces
coincident with refinement boundaries without assigning a film to the solid side.
The original nodal lattice may be an integer isotropic refinement of the current
finest tree, so it can be retained after coarsening. EB option presence, lattice
dimensions and every nodal value must agree across ranks; empty ranks participate.
Owner/capacity checks and all payload exchanges remain collective and transactional.
Probe endpoints must stay inside the root domain, with the same offset and
normal/location tolerances as serial EB film exchange. Films have no global ID;
receive ordering remains source rank then local order. This adds migration only;
distributed EB impact, film exchange and checkpoint integration remain open.
See [MPI film ownership evidence](validation/spray-mpi-eb-film-20260918.md).

`advance_mpi_amr_film_coupling` now accepts optional `eb_geometry` and
`nodal_phi` together. It uses EB film migration, fluid-side sampling and the
volume-weighted composite EB gas depositor with common adaptive film substeps.
The whole interval retains transactional owner-local gas, temperature, films,
component inventories and global wall heat. EB mode requires nonperiodic outer
faces and one spatially uniform wall temperature (all six entries of the
existing `spray_wall_options%temperature` must match). The existing wall option
validator is applied with slip-wall face types; enabled flags must satisfy its
evolution contract but do not select EB facets. Material-specific conductivity
models remain available. Geometry/EB option presence is agreed collectively;
model/property tables remain identical across ranks by caller contract.
Gas feedback preserves covered cells and restricts hidden parent fluid state.
This source API does not yet connect distributed particle impacts, reactive
composition or the MPI application/checkpoint path. See
[MPI EB film exchange evidence](validation/spray-mpi-eb-film-step-20260918.md).

`advance_mpi_eb_spray_reactive` now accepts optional `films`, `film_mass` and
`film_settings` together, plus an optional `film_steps` output. It composes
particle half / film half / reactive gas / film half / particle half intervals.
All optional-state presence is agreed across ranks before entering collectives.
Gas, particles, films, component arrays and diagnostics publish only after the
entire composition succeeds; wall heat joins the existing gas and particle wall
ledger. `particle_steps` continues to count only particle substeps, while
`film_steps` counts both film halves and is zero on failure. Film settings retain
the uniform-temperature and nonperiodic EB source restrictions above. Without
`last_id`, parcels retain specular reflection and films come from the caller.
Supplying the optional shared `last_id` with all three film arguments enables
EB impacts in both particle halves. The watermark is private until the entire
step succeeds, together with new children and films; optional `model_rejected`
propagates the source diagnostic. Films created in the first particle half can
evolve in the following film halves. Films created in the final particle half
begin thermal evolution on a subsequent step. Presence of the watermark and
film arguments must agree across ranks. The MPI application connects this path
as described below. See [reactive MPI film evidence](validation/spray-mpi-eb-film-reactive-20260918.md)
and [reactive impact evidence](validation/spray-mpi-eb-impact-reactive-20260918.md).
A watermark may also be supplied without films for KH children with specular
flight. KH mode requires it; its presence is agreed across ranks and both
particle halves use a private copy. It commits only after the complete reactive
step succeeds. The MPI executable now enables this path and restores injection
plus four-child ID groups from its existing checkpoint payload. See
[KH application evidence](validation/spray-mpi-eb-khrt-application-20260918.md).

The MPI EB reactive wrapper also applies `options%maximum_multiplicity` after
the final particle half, before publication. A positive cap requires shared IDs
and must agree across ranks. Splitting failure preserves the entire incoming
state; optional `split_created` reports newly allocated IDs only on success.
The application saves their cumulative count for strict restart ID accounting.
See [computational splitting evidence](validation/spray-mpi-eb-split-20260918.md).


`spray_mpi_eb_impact_mod` exports `impact_mpi_level_set` for one earliest contact
per supplied trial trajectory. It uses original nodal geometry, sampled parcel
pressures and common EB wall settings, queries pre-event film thickness across
ranks without collecting the films, and calls the existing collective impact
resolver/ID allocator. Child IDs use rank-prefix reservation above the shared
watermark. `hit` and `fractions` refer to the original population; the resulting
population includes appended children and retained zero-mass parents. `wall_change`
is the globally summed inventory received by the wall. All payloads, IDs and
diagnostics are unchanged on failure except `ok`/`model_rejected`. Invalid
empirical output is reported collectively before publication.
Bounds, nodal dimensions and every nodal value must agree across ranks;
model/property tables are identical by caller contract. This is a batched
contact API: newly deposited films do not affect other contacts in the same
batch. The caller owns remaining flight, dead-parent removal, final particle/film
migration and binding the nodal field to the gas geometry. The source and reactive integrations are described below;
application integration is described below. See
[MPI EB contact evidence](validation/spray-mpi-eb-impact-20260918.md).

`advance_mpi_level_set_impacts` now wraps collective contacts in a force-free
remaining-flight interval. It inherits frozen pressure, remaining time and
contact depth from each parent to its children, repeats collective batches
until every live trajectory finishes, and compacts consumed parents only after
success. Empty/inactive ranks participate until global completion. ID allocation,
film patches, parcel/component arrays, global wall receipt and global hit count
publish together; any failure retains their input values. The contact limit
applies along parent/child trajectories, including consumed intermediate parents.
Film thickness is frozen within a batch and updated for subsequent batches;
contacts are not globally sorted by physical time. Parcels/films retain their
source-rank ownership during this interval. Gas feedback, film evolution and
final finest-owner migration are the caller's responsibility. See
[remaining-flight evidence](validation/spray-mpi-eb-impact-interval-20260918.md).

Its optional `outflow_change(5+nc)` enables all-root-face outflow. Presence must
agree across ranks. Each flight is clipped at its first exterior crossing;
earlier EB contacts use the corrected remaining time. Escaped records are
removed without recycling IDs. The global negative system inventory contains
mass, momentum, energy and component masses in component order, while the wall
output remains wall receipt. Exterior events share the per-lineage hit budget.
All ranks receive the same inventory only after collective success. The tree
migration wrapper `advance_mpi_eb_tree_impacts` accepts the same optional output
and retains it privately until final particle and film migration both succeed.
MPI source, reactive coupling and the MPI application now expose the option below. See [MPI outflow evidence](validation/spray-mpi-eb-impact-outflow-20260918.md)
and [migration evidence](validation/spray-mpi-eb-outflow-migration-20260918.md).

The MPI AMR source API accepts `eb_outflow_change(nvar)` with nodal EB geometry,
using either specular reflection or impact films. Presence enables all-face particle outflow and must agree
across ranks. Component-order losses are mapped into vapor-species entries;
already-global impact inventories are accumulated without a second MPI sum.
The existing `boundary_change` continues to report EB system wall exchange.
Both ledgers publish only after the complete source interval and final ownership
migration succeed. KH eligibility clips the candidate flight at the root before
querying the EB. With specular reflection, local drift losses are summed once
across ranks before accumulation, and escaped records are removed after each
accepted substep. This prevents later gas sampling at their exterior endpoints.
The reactive wrapper and application support both specular and impact-film outflow. See
[specular source evidence](validation/spray-mpi-eb-specular-outflow-20260918.md) and
[source evidence](validation/spray-mpi-eb-outflow-source-20260918.md).

`advance_mpi_eb_spray_reactive` accepts `particle_outflow=.true.` with specular
reflection or with impact films and the shared ID watermark. The effective flag must agree across
ranks. Both particle halves collect their exterior inventories privately; these
join the gas boundary output only after all particle, film, reaction, transport,
hydro and optional splitting work succeeds. Wall exchange remains separate.
See
[reactive evidence](validation/spray-mpi-eb-outflow-reactive-20260918.md).
The MPI EB spray application accepts the same `particle_outflow` flag in
`eb_spray`, with either wall model. It enters the serialized restart
context. Existing particle/film and gas archives preserve surviving sparse IDs,
watermark and boundary history through rank-changing restart. See
[application evidence](validation/spray-mpi-eb-outflow-application-20260918.md).





`advance_mpi_eb_tree_impacts` composes that interval with transactional AMR
ownership. It first validates/migrates a private copy of existing EB films,
advances the impact interval with pressures still attached to their input
parcels, then migrates live parcels and films to the finest current fluid-side
owners. Failure at either final receive-capacity check restores the entire
pre-impact state, including IDs and diagnostics. The root upper bound is derived
from tree shape and root spacing. Both capacities come from the supplied particle
and wall options. This wrapper adds final migration, not gas feedback or film
evolution; callers must resample pressures after repartitioning for the next
source interval. See
[impact/migration evidence](validation/spray-mpi-eb-impact-migration-20260918.md).

`advance_mpi_amr_mixture_coupling` also accepts its existing film arguments and
shared `last_id` with `eb_geometry`/`nodal_phi`. It advances gas/parcel drag,
heating and evaporation, then collective EB impact and remaining child flight
on each accepted source substep. Existing films are validated/migrated initially;
live parcels and deposited films migrate to final finest fluid-side owners.
The global impact receipt enters `boundary_change` once with the system-gain
sign. Gas, temperature, particles, films, both component matrices, IDs and the
ledger publish together after the full interval succeeds. Optional
`model_rejected` reports inadmissible empirical output collectively without
repeated timestep halving.

EB impact mode requires one uniform wall temperature across the six existing
wall-option entries; enabled face flags do not select EB facets. Contact depth
is limited to 32 per source substep. KH child shedding is now supported by the
MPI source API, with wall-crossing precedence as described below. Properties/settings must match across ranks. Film thermal evolution
is still supplied by the separate film source API; these source arguments alone
do not evolve films. The MPI reactive wrapper can compose these impact-enabled
particle halves with separate film evolution, as described above. MPI application
and restart controls are described below.
See [source coupling evidence](validation/spray-mpi-eb-impact-source-20260918.md).


`spray_eb_field_mod` supplies `eb_spray_stencil`, `sample_eb_spray_field` and
`deposit_eb_spray_gas`. Callers must supply geometry and its original nodal level
set together. Point queries reject solid locations. Up to eight neighboring
cell-center weights are multiplied by fluid volume fraction, restricted to cells
connected to the containing cell by open faces within that support, then
normalized. The same weights sample fields and deposit integrated contributions;
deposition divides by each cell's actual fluid volume. This is a bounded
volume-weighted cloud-in-cell approximation, not exact cloud/polyhedron overlap
or a linear-exact EB interpolant. Exterior support is truncated and normalized.
The gas depositor validates incoming active states and recovers final EOS
temperatures before publishing state and temperature together. Covered cells
are untouched. Sampling also publishes only after every query succeeds. The
caller must provide consistent mass/species feedback; this layer has no small-cell
redistribution or timestep selection.

`spray_eb_step_mod`'s `advance_eb_mixture_coupling` connects these APIs to the
existing multicomponent drag/heating/evaporation kernel and representative-drop
ETAB/KH-RT breakup. It selects common particle substeps from drag, component
loss, temperature change and cell transit limits. Each substep applies sources
first, then follows the updated velocity through stationary specular EB contacts;
this is a first-order source/drift split. The source kernel's free-space
trapezoidal position update is replaced by this EB trajectory. Feedback is
deposited at the pre-substep position. Wall impulse/energy is returned separately
as `boundary_change(5)` with the particle-change sign; it is not deposited into
gas. Failed gas recovery or trajectories retry up to 30 halvings, and the entire
interval publishes gas, temperature, parcels, component masses and wall ledger
together. Failed intervals return zero substeps and preserve other outputs.
The API retains zero-mass records, bounds the population to 125000 and rejects
KH child creation. Each drift allows at most 32 contacts and requires trial
endpoints within the root domain. Exterior boundary events, liquid wall films
and MPI/AMR evolution remain separate work.

`spray_eb_reactive_mod` adds `advance_eb_spray_reactive`: a particle half-interval,
the existing `advance_reactive_eb_full_3d` gas update, then a second particle
half-interval. It passes gas solver, redistribution, chemistry integrator and
molecular-transport controls to the existing EB operator. Particle exchange
remains first order; the symmetric outer composition alone does not establish
second-order accuracy. All fields, component masses and diagnostics publish
together only after both particle halves and the gas step succeed. Returned
ledgers separate exterior boundary exchange, combined gas/particle EB wall
exchange, and gas chemistry; each uses the change-in-system sign convention.
`particle_steps` counts both accepted particle halves and is zero on failure.
The allocated interface additionally applies computational splitting after the
last half, with a shared ID watermark and optional `split_created` output. The
fixed-population interface rejects an enabled multiplicity cap. The serial
application uses the allocated interface on both wall paths and persists cumulative
split IDs for restart accounting. The shared `split_spray_parcels` helper accepts
optional `watermark` and `component_mass` arguments and publishes them atomically
with descendants. See [serial splitting evidence](validation/spray-eb-split-20260918.md).
The caller must select a stable gas timestep; this wrapper does not retry or
adapt gas timesteps. An experimental fixed/selected
[serial EB spray application](../cases/eb_spray_3d/README.md) now connects nodal
geometry, scheduled cone injection, coupled stepping with bounded retries, gas and
particle CSVs and cumulative boundary/wall/chemistry ledgers. Immutable serial
checkpoints now preserve gas, particle/component state, event/ID counters and all
cumulative ledgers; restart checks the complete numerical and nodal geometry
context. Up to 1024 ordered injection events share
up to 16 cone ports with per-port geometry, population and composition, per-event
mass scaling and zero-dose support. Nozzle models, size distributions, rate
profiles, outlet radii and discharge coefficients accept per-port overrides. Steps
end at event times; tied and final-time events are included. Fluid-location
validation precedes population/ID/injected-inventory publication. See
[serial port evidence](validation/spray-eb-ports-20260918.md).

The MPI tree module now exports `migrate_eb_tree_spray_parcels`, a nonperiodic
ownership primitive. It accepts regular tree topology/owners and the replicated
root-domain nodal level set at the finest tree resolution. The field must be the
same one used to build the gas EB hierarchy. Private candidates pass through
existing AMR migration, then the routine collectively checks nodal dimensions,
exact rank-to-rank field agreement and fluid point locations before publishing
parcels/components. A solid point or any other rank's failure preserves all
caller populations. It retains the regular half-open outer-domain ownership
rule. This is nodal point classification, not validation of the gas hierarchy's
rounded/coarsened volume metrics. It does not yet provide distributed EB
sampling, deposition, time integration or restart. The replicated geometry is
bounded to one million finest cells; full-field agreement is not a scalable
geometry distribution strategy.

`spray_eb_tree_stencil_mod` now provides `eb_tree_spray_stencil` for composite
AMR support. It takes topology, one geometry per node and the finest root-domain
nodal field. Existing rectangular-cloud overlap supplies exposed coarse/fine
cells; their weights and physical volumes are multiplied by fluid volume
fraction. A segment query from the particle to each fluid centroid removes
support hidden behind solid geometry. Remaining weights are normalized and
zero-weight entries removed. Solid particles, inconsistent physical patch bounds
or an empty visible support fail without changing caller arrays. This is a
volume-fraction-weighted, centroid-visibility approximation, not exact cloud/EB
intersection or a linear-exact coarse/fine interpolant. The field and node
metrics must come from the same EB hierarchy. The stencil alone does not sample,
deposit, redistribute small-cell sources or perform MPI field exchange.

`sample_mpi_amr_tree_spray_field` accepts paired optional `eb_geometry` and
`nodal_phi` arguments to use these EB stencils in the existing owner-routed
request exchange. Both must be present on every rank and all axes nonperiodic.
Replicated nodal values, node bounds, fluid fractions and fluid centroids are
checked for exact rank agreement, including ranks without queries. Gas fields
remain owner-local; only requested cell values travel. All local sample arrays
publish after collective success. Request-capacity, geometry and solid-point
failures preserve existing samples.

The field and gas deposition APIs now accept the same paired EB arguments.
Owner-routed contributions divide by physical fluid volume. Parent restriction
communicates averaged `volume_fraction * state` and volume fraction, then
divides to obtain fluid-volume-weighted parent values. Only parent cells hidden
by children are replaced; exposed and solid cells retain their deposited/original
state. Restricted fractions must agree with parent EB metrics to roundoff.
The gas wrapper recovers EOS temperatures only in positive-fluid-volume cells
and commits gas and temperature together after all ranks succeed. Failed source
updates, inconsistent coarse/fine fluid volumes or invalid EOS states preserve
the whole incoming gas/temperature state. Small-cell source redistribution is
not provided; timestep rejection remains the caller's responsibility.

`advance_mpi_amr_mixture_coupling` accepts paired `eb_geometry` and `nodal_phi`
to connect MPI EB sampling/deposition, mixed-liquid drag/heating/evaporation,
representative-drop breakup, specular drift and final finest-owner migration.
Common source/cell-transit substeps and bounded trial halving remain active.
The EB path replaces the free-space position update with source-first drift at
updated velocity, as in the serial EB interval. Local wall exchanges accumulate
only for accepted substeps and are summed once into the returned global boundary
ledger. The full interval publishes gas, temperature, particles, component masses
and ledger together. Exhaustion or late migration failure leaves them unchanged.
EB mode rejects explicit exterior-face settings; trajectories must stay inside
the root domain. Optional wall-film state enables the impact path described above.
KH child creation requires a shared `last_id`. Each parent trajectory is checked
against the original EB field before release: a crossing suppresses that parent's
KH update for the source substep, leaving the wall path to handle it. Unhit
parents can release four children through the collective allocator, with a gas
energy debit for their transverse motion. Parents and new children then undergo
the selected specular or impact remaining-flight path and final owner migration.
Source-first breakup/flight splitting remains in force; this is not a globally
time-ordered aerodynamic-breakup event solver. The collective KH helper's
optional `eligible` mask preserves excluded parent mass, pending breakup state
and zero release-energy contribution. See
[KH source evidence](validation/spray-mpi-eb-khrt-source-20260918.md).
The MPI application now connects KH mode and its shared IDs as described below. This is a particle
interval; the distributed gas/spray wrapper below composes it with gas dynamics.

`spray_mpi_eb_reactive_mod` exports `advance_mpi_eb_spray_reactive`: particle
half-step, existing MPI EB reaction/transport/hydro step, particle half-step.
It accepts the same owner-local gas/temperature tree and replicated geometry
and original nodal field. Gas and particle state, component masses, all three
global ledgers and the transport limiter publish only after all stages succeed;
failure returns zero particle steps and preserves the caller state/diagnostics.
The wall ledger includes both particle impulses; exterior and chemistry ledgers
come from the gas operator. All use change-in-system signs. Optional gas
integrator and subcycling settings pass through to the existing gas solver.
The particle restrictions above remain. Fixed and selected MPI EB spray
executables now connect scheduled mixed-liquid cone injection, coupled evolution,
combined gas/particle budgets and per-rank particle CSV output. See the
[example and controls](../cases/mpi_eb_spray_3d/README.md) and
[application evidence](validation/spray-mpi-eb-application-20260918.md).
Scheduled events clip the common step to the next injection time; tied events,
zero scales and final-time events are supported, with global live-population
limits and monotonic IDs. The budget has a separate injected column and starts
from gas before the first event. See
[schedule evidence](validation/spray-mpi-eb-schedule-20260918.md).
A rank-neutral particle payload API, `mpi_eb_particle_checkpoint` in
`spray_mpi_eb_particle_checkpoint_mod`, now preserves parcels, liquid component
masses, time/step, next event, ID watermark and injected inventory. It gathers
private copies for immutable publication and restores onto the current finest
fluid owners transactionally. The caller must bind its single-line context to
the gas checkpoint generation, geometry, model and schedule; this payload does
not itself contain or verify gas state. See
[payload evidence](validation/spray-mpi-eb-particle-checkpoint-20260918.md).
Its optional `wall_films`, `film_mass`, `maximum_films` group now stores stationary
EB films in `PELEF_EB_PARTICLES_V2`. Calls without that group retain the original
V1 format. Readers require the format corresponding to their argument group;
film-mode presence must agree across ranks. Film gathering reuses the physical
EB tree with every node temporarily assigned to the archive writer; the film
capacity therefore also bounds the gathered archive. Reads parse and validate
private film payloads, then migrate them to current fluid-side owners before
publishing particles, films or metadata. Invalid film geometry after particle
redistribution rolls back the whole restore. Film ordering follows migration
source order; films have no global IDs. See
[film payload evidence](validation/spray-mpi-eb-film-checkpoint-20260918.md).
The MPI application now selects V2 when wall impacts are enabled.

The fixed/selected MPI EB spray applications now connect this payload with the
existing gas snapshot using a generation-bound manifest. Restore stages both
states and verifies the complete coupled budget and injection position before
publication. Immutable step-suffixed checkpoints support changed-rank continuation;
see the [case guide](../cases/mpi_eb_spray_3d/README.md) and
[restart evidence](validation/spray-mpi-eb-restart-20260918.md).
Whole-patch rebalancing now uses the existing measured gas work and moves
owner-local gas/temperature and particles atomically. `rebalance_mpi_amr_spray`
accepts optional `nodal_phi` for EB owner migration; this mode requires
nonperiodic coordinates and accepts the existing `wall_films`, `film_mass`,
`maximum_films` group. EB films migrate using fluid-side probes against the
original nodal field. A rejected particle or film move preserves the proposed
gas/owner changes privately. See
[rebalancing evidence](validation/spray-mpi-eb-balance-20260918.md).
`remap_mpi_eb_spray` in `spray_mpi_eb_remap_mod` now stages the existing
conservative MPI EB gas remap and particle owner migration before publishing
topology, owners, geometry, gas/temperature, parcels and liquid components.
Movement of refined patches, coarsening and refinement share the same physical
root domain and immutable geometry. The original nodal field may be a common
integer refinement of the current finest tree lattice, allowing coarsening
without replacing the original collision geometry. MPI sampling requires exact
nodal dimensions and values across ranks. See
[remap evidence](validation/spray-mpi-eb-remap-20260918.md).
The optional `wall_films`, `film_mass`, `maximum_films` group also carries EB
films through remap. Presence must agree across ranks. Film migration occurs on
private state after gas and particle remap; any geometry or receive-capacity
failure discards the entire proposed topology/state change. Film locations,
normals and physical/component payloads remain unchanged. Both remap and
rebalancing retain their particle-only calling forms. The library connection is
qualified separately from MPI application film persistence; see
[film repartition evidence](validation/spray-mpi-eb-film-repartition-20260918.md).
The fixed/selected MPI EB spray applications now route prescribed and
existing temperature-tagged regrids through this coupled remap. Checkpoint
restore prepares the stored hierarchy privately, checks the saved stage/time,
then publishes hierarchy, gas and particles together. Original nodal values
remain fixed across coarsening/refinement. See
[application regrid evidence](validation/spray-mpi-eb-regrid-application-20260918.md)
for prescribed schedules; adaptive tagging is not newly qualified there.

The MPI EB spray application now accepts up to 16 independent cone ports and an
`event_ports` selector in the shared ordered event table. Ports may differ in
count, base mass, position, direction, angle, diameter, temperature and liquid
composition. Omitted port entries inherit the existing single-port settings.
Component property tables remain shared; per-particle component masses determine
subsequent mixture evolution. The particle material ID records its injection
port, including across remap, migration and restart. See the
[port input guide](../cases/mpi_eb_spray_3d/README.md) and
[port evidence](validation/spray-mpi-eb-ports-20260918.md).
Port-specific 2--64 point mass-flow tables now use the existing piecewise-linear
integrator over explicit event start/end intervals. Their integrated masses,
scaled by the event multipliers, replace prescribed packet masses. Zero-area
events reserve no IDs; restart reconstructs the effective positive-event counts.
The complete profile/interval input is context-bound. See
[profile evidence](validation/spray-mpi-eb-profile-20260918.md). See
[coupled-step evidence](validation/spray-mpi-eb-reactive-20260918.md).

The `pelef_spray_tree_runtime` library (or
`pelef_selected_spray_tree_runtime` for a configured bundle) exports
`locate_amr_tree_spray_cells` and `partition_amr_tree_spray_parcels` from
`spray_amr_tree_partition_mod`. Inputs use the physical root lower corner and
root cell spacing, a periodic flag per axis and zero-based node owners. Lookup
returns `(node,i,j,k)` with one-based local cell indices. Cells and patches
include their lower faces and exclude upper faces; periodic coordinates wrap
for lookup without changing parcel positions. Nonperiodic exterior points fail.
The partitioner requires strictly increasing input IDs and preserves that order
within each destination rank, returning counts and zero-based offsets. Its
location array follows the reordered parcels. All outputs retain their previous
values on validation failure. Empty populations/ranks are supported. Coordinates
follow the existing spray bound of magnitude at most 1e9 m; root spacings are
between 1e-12 and 1e9 m, and integer level domains are limited to 2^52 cells per
axis for floating coordinate lookup. This API does not classify EB solids.

With MPI enabled, the tree library exports `migrate_amr_tree_spray_parcels`
from `spray_mpi_amr_tree_mod`. Every rank supplies the same validated topology,
owner map, root bounds/spacing and periodic flags; component payload presence
and count must agree. Input parcel IDs are locally increasing and must be
globally unique. Received records are sorted by ID, and duplicate IDs arriving
at one rank reject the transaction. Each rank can set its own receive capacity
up to one million parcels. Capacity, metadata or received-state validation
failure leaves all ranks' original parcel/component arrays unchanged. Empty
populations/ranks, periodic lookup, ownership changes and regrid reassignment
are supported. This API migrates an already accepted particle population; it
does not update gas fields, particle trajectories, checkpoint data or EB impacts.

`spray_amr_tree_coupling_mod` supplies `amr_tree_spray_stencil`,
`sample_amr_tree_spray_field` and `deposit_amr_tree_spray_field`. A rectangular
particle cloud has the width of its finest containing cell in each axis. Its
overlap with exposed composite cells gives positive normalized weights; covered
parent cells receive no direct weight. On a uniform mesh this matches the
existing trilinear cloud-in-cell interpolation. At coarse/fine interfaces the
field is treated as piecewise constant inside each exposed cell; this is not a
linear-exact coarse/fine reconstruction. Periodic cloud pieces wrap, while
nonperiodic clouds clip to the domain and renormalize. At most one million
overlap entries are accepted (the stencil API permits a smaller capacity).

Sampling and deposition use the same weights. Deposition accepts integrated
amounts, divides by each receiver's physical cell volume, and averages covered
parents after all leaf updates succeed. Outputs commit only after success;
empty deposits preserve the entire input field. The generic field API does not
enforce gas thermodynamics. Its `deposit_amr_tree_spray_gas` wrapper validates the
initial gas, applies the same conservative field update, then recovers EOS
temperatures on every patch, including averaged parents. Conserved fields and
temperature fields commit together only if every cell accepts the update.
Invalid density, species composition or internal energy reject the entire
feedback. Empty gas feedback validates inputs and preserves both fields. This
is the serial tree API; its owner-local MPI counterpart is described below.
The MPI particle interval described below uses this feedback design. No AMR
convergence or physical spray validation is claimed.

`sample_mpi_amr_tree_spray_field` in `spray_mpi_amr_field_mod` applies the same
stencil to owner-local fields. The caller supplies the common component count,
tree, owner map, physical bounds/spacing and periodic flags. Particle positions
are local to each rank and can be empty even on ranks serving gas queries.
Only requested cell values are fetched from owners, with no full-field gather;
covered/unrequested payloads are not scanned. Nonfinite requested values,
invalid local positions, field-shape/metadata disagreement or request-capacity
failure reject collectively and preserve every rank's previous sample array.
The optional `maximum_requests` bounds each rank's total cloud-cell requests
(default one million, range zero to one million); ranks may choose different
capacities. Received values are accumulated in the serial stencil order.

`deposit_mpi_amr_tree_spray_field` uses the same distributed stencil preparation
and routes weighted, volume-scaled contributions directly to cell owners.
Each owner applies source ranks in increasing order, retaining each source's
particle/stencil order, then the existing MPI restriction refreshes covered
parents. No full-field gather is used. All owned field values must be finite.
The optional `maximum_requests` bounds outgoing cloud-cell entries and
`maximum_received` independently bounds incoming entries on each rank (both
default to one million, with range zero to one million). Invalid contributions,
capacity or arithmetic failure preserve the original distributed field.
Globally empty deposition validates inputs and leaves the field unchanged.
For a fixed source ordering, deposition matches serial accumulation; changing
particle ownership can change summation roundoff. This is a generic density-field
operation without thermodynamic constraints. Its `deposit_mpi_amr_tree_spray_gas`
wrapper validates the initial gas, deposits into private owner-local fields and
recovers temperatures collectively after parent restriction. U/T publish only
when every owner accepts every cell. A failure on any rank preserves all ranks'
original U/T, and globally empty feedback validates without changing either
field. The caller supplies an identical ordered species table on every rank;
the wrapper retains the request/receive-capacity controls of field deposition.
`advance_mpi_amr_mixture_coupling` in `spray_mpi_amr_step_mod` connects
mixed-liquid particle intervals. It samples conserved fields and temperature
guesses together, then uses the existing drag/heating/evaporation and breakup
kernels. Every rank takes a common minimum substep; the finest spacing anywhere
in the tree bounds particle CFL, including approaches to refined regions.
Failed particle updates or thermodynamic feedback halve the trial interval, up
to 30 attempts. ETAB/KH-RT representative breakup and optional KH children use
the existing component inventories, collective IDs and gas-energy debit.
Feedback is deposited at the pre-step particle positions. By default accepted
trajectories wrap periodically; final survivors migrate to their finest-patch owners.
Gas, temperatures, parcels, component masses and optional ID watermark publish
only after the entire interval and final migration succeed. Substep exhaustion
or receive-capacity failure therefore rolls back the whole interval.

All ranks must supply the same species, transport, liquid tables and options.
Particle IDs must be globally unique and locally increasing. The interval and
component count are checked collectively; local particle populations may be
empty. This API advances particle exchange against the gas, without gas
hydrodynamics/chemistry, regridding, injection or checkpoint orchestration.
Planar wall-film coupling is described below; EB collisions remain open.

The particle-interval API also accepts optional `boundary_faces` and
`boundary_change`. Periodic faces wrap; inflow/outflow faces remove exiting
parcels, and slip/no-slip/isothermal wall faces use the existing specular parcel
reflection model. Gas wall names here do not imply liquid deposition, thermal
wall exchange or film formation. Face types must agree across ranks. The global
integrated boundary ledger contains gas-plus-particle inventory changes in
mass, momentum, total energy and species order. It includes represented liquid
mass/momentum/energy carried out of the domain and wall impulse. The ledger
publishes only with whole-interval success and otherwise retains its incoming
value. The reacting wrapper and application accept the same faces and ledger;
checkpoint continuation preserves the cumulative inventory.

The particle interval accepts optional `wall_films`, `film_mass` and
`film_settings` together, requiring explicit faces and a shared `last_id`.
It reuses `impact_mpi_planar_walls` before generic boundary handling: earliest
nonperiodic crossing, distributed existing-film thickness, mixture impact
classification, deposited films and optional splash children. KH children use
the parent's pre-step position and sampled pressure for impact geometry.
Accepted substeps retain private films and accumulate the already global wall
impulse/energy separately from rank-local flow exchanges. Final migration and
publication include films; the returned boundary ledger subtracts what the wall
received. Rejected trials and whole intervals preserve incoming film/ID state.
This retains the existing planar impact model, not EB collision geometry.

`advance_mpi_amr_spray_reactive` in `spray_mpi_amr_reactive_mod` composes a
particle half-interval, the existing AMR reaction/transport/hydro split, and a
second particle half-interval. The gas step accepts its existing reconstruction,
limiter, Riemann solver, chemistry integrator, transport switches and optional
recursive subcycling. Optional `boundary_faces` defaults to periodic;
`inflow_state`, `inflow_temperature`, `wall_temperature` and `wall_velocity`
forward the existing gas boundary data. Optional `boundary_change` returns the
sum of both particle halves and the gas boundary inventory, replacing its input
only on success. It uses conserved mass, momentum, total energy and species
order; chemical species conversion is not a boundary contribution. The wrapper keeps all stage
outputs private: gas, temperatures, particles, component masses, ID watermark
and gas correction/limiter diagnostics publish only after both particle halves
and the gas step succeed. Failure returns zero particle substeps and leaves
caller state and diagnostics unchanged. The caller chooses/retries the overall
gas timestep; particle halves select their own common substeps.

Optional `wall_films`, `film_mass` and `film_settings` must be supplied together
and consistently on all ranks, with `last_id` for impact-created children.
They wrap the composition in film half-intervals:
F/2–S/2–gas–S/2–F/2. Film gas exchange and wall heat use the existing AMR film
interval, with private film state until the full step succeeds. The boundary
ledger includes both film wall-energy contributions. Optional `film_steps`
reports only accepted film substeps; `particle_steps` retains its particle-only
meaning. Failure returns zero counts and preserves films along with the other
state and diagnostics. Existing calls without films retain the original split.
Particle halves also create films through the existing planar impact model.
The fixed/selected AMR applications expose this through `walls` settings and
carry films through adaptation, rebalancing and checkpoint continuation.

The MPI tree spray library now links the existing MPI tree gas runtime instead
of compiling separate copies of its cell-query/restriction/thermodynamic modules.
The selected library links the selected MPI tree runtime; it does not acquire
the fixed-mechanism core. The MPI AMR spray application below supplies setup,
initial or timed injection, prescribed or automatic regridding and checkpoint/restart.

`remap_mpi_amr_spray` in `spray_mpi_amr_remap_mod` connects conservative gas
remapping and particle migration for a changed regular tree or owner map. It
validates the original gas thermodynamics, migrates private particle/component
arrays to the target finest owners, then remaps gas and recovers temperatures.
Tree, owner map, gas, temperatures, parcels and component inventories publish
together only after every stage succeeds. The root resolution and physical
domain must remain fixed; refinements may move, appear or disappear. Particle
positions and breakup/lineage payloads remain unchanged. Periodic flags control
lookup; nonperiodic exterior parcels reject. The caller retains the unchanged
ID watermark and chooses the new topology/owners and per-rank receive capacity.
This API does not choose the new topology. The MPI AMR spray application can
schedule one caller-prescribed topology file. Optional `wall_films`, `film_mass`
and `maximum_films` must be supplied together and consistently on all ranks.
Their private finest-owner migration follows successful gas/particle remap, so
film capacity or geometry failure also preserves the original gas, particle,
film and topology state. Exact wall coordinates and all film fields survive
coarsening/refinement. Film-free calls retain their existing behavior.

`pelef_mpi_amr_spray_3d` and `pelef_mpi_amr_spray_3d_selected` share
`spray_mpi_amr_application_mod`. They take gas, topology and spray input paths,
initialize owner-local gas, inject mixed-liquid cones from up to 16 ports at specified times with collective IDs,
advance coupled steps with whole-step retries and remap at a prescribed time
or at step-based temperature-adaptation checks. Gas CSV contains only exposed cells; separate per-rank particle
CSVs are analysis outputs. See the [case guide](../cases/mpi_amr_spray_3d/README.md)
for controls and limitations. Gas-work-based rebalancing moves owner-local gas
and particles together after accepted steps. Gas per-face reservoir inflow,
outflow and adiabatic/isothermal tangentially moving walls now connect to the
particle removal/specular model. Regridding, tagging, balancing and injection
migration use the corresponding periodic axes. Reservoir CFL and transport
bounds use the current finest spacing after adaptation. Optional `walls` controls
enable planar impacts and film evolution, including film CSVs and restart; EB
particle handling remains open.

`migrate_mpi_amr_wall_films` in `spray_mpi_amr_film_mod` assigns stationary
mixed-liquid film patches to the owner of their finest containing AMR cell.
The physical domain is closed for this lookup: exact upper-wall coordinates
select the last cell at each level, while stored positions/normals remain exact.
Exterior patches reject rather than wrapping. The API preserves mass,
temperature, diameter, height, material and every component mass. Films have no
global ID; receive order is source rank followed by original local order.
It checks collective topology/owners/domain/component count/capacity and stages
all output until validation and exchange succeed. The one-million-per-rank
send/receive bound matches the existing Cartesian film migration. Empty ranks
and globally empty films are supported. This is a library ownership operation;
Particle intervals create films; the application passes them through remap,
balancing and checkpoint operations when wall impact is enabled.

`rebalance_mpi_amr_spray` also accepts the optional group `wall_films`,
`film_mass`, `maximum_films`. When gas-work costs select new owners, films
migrate privately after gas and particles. All inventories and owners publish
together; failed film migration preserves the original state and returns
`changed=.false.`. Film arrays/settings presence agrees across ranks and film
inventories are checked even when costs select no ownership change. A no-change
result leaves film ordering and storage unchanged. Costs still represent gas
work; this does not add film/particle-aware cost estimation.

`advance_mpi_amr_film_coupling` in `spray_mpi_amr_film_step_mod` advances
stationary planar mixed films with existing heating/evaporation kernels and
conservative owner-local gas feedback. Film normals must select an enabled
axis-aligned physical wall and their normal coordinate must equal that wall.
Closed-wall positions stay exact; gas probes clip to a half finest-cell spacing
inside the domain before using the existing composite stencil. Sampling combines
conserved state and temperature guesses; deposition averages covered parents
and recovers EOS temperatures. All ranks use a common minimum substep from
component depletion and temperature change, with bounded halving on failure.
Gas, film/component arrays and the integrated wall-energy ledger publish only
on whole-interval success. Failed intervals preserve the incoming ledger and
state. This advances exchange against gas, not gas transport/chemistry or impact
creation; the reacting wrapper supplies coupled splitting. Application film
controls are described in the AMR case guide.

`mpi_amr_spray_checkpoint` in `spray_mpi_amr_checkpoint_mod` writes distributed
indexed gas payloads and a rank-zero particle archive, then publishes an immutable
completion marker with a POSIX hard link. Both payloads bind the same generation
and numerical/model context. Reads validate both before publishing topology,
owner-local gas and redistributed particles. The archive includes composition,
lineage and breakup fields, ID watermark, time/steps, next injection event,
regrid status, injected inventory and cumulative boundary inventory. The latter
uses the existing gas checkpoint budget field without changing particle records.
Optional checkpoint `periodic` axes bind nonperiodic lookup to the context;
optional `boundary` publishes only after successful read. Omitting this output
rejects a saved nonzero budget instead of silently losing it. Rank-zero particle collection is bounded
to one million global parcels. The application exposes one checkpoint time,
optional stop and changed-rank continuation; see the case guide for file handling.

Optional checkpoint `wall_films`, `film_mass`, `maximum_films` arguments form
one collectively agreed group. In film mode the particle payload uses
`PELEF_AMR_PARTICLES_V2`, with a bounded film count and explicit material,
position, normal, mass, temperature, diameter, height and component records after
the parcel records. The numerical/generation context binds film presence. The
manifest format and file set stay unchanged. Film-free payloads remain V1.
Films gather privately on rank zero, with a one-million-global-film bound.
Reads validate film records and the final footer/EOF, then migrate to finest
owners on the reading communicator. Gas, particles, films and all metadata
publish together only after every validation and capacity check succeeds.
A different rank count, including one rank, is supported by the library. The
application passes these arrays when any `walls%enabled` face is active.


## PPM reconstruction

The fixed/selected serial spray application accepts
`hydro_reconstruction='characteristic_ppm'` in `&spray`. PPM uses the existing
particle/gas split, and checkpoint context binds the reconstruction selector.
A fixed Mach-inlet fixture with molecular transport, Smagorinsky and evaporating
parcels has exact restart evidence
([focused checks](validation/ppm-regular-spray-20260925.md)). The regular MPI spray application also accepts this selector in `&mpi_spray`,
with fixed changed-rank continuation
[checks](validation/ppm-regular-mpi-application-20260925.md).
MPI AMR spray instead uses the gas namelist's `reconstruction` selector.

For serial fixed/selected spray, `difmag=-1` (the default) selects `0.1`
for all supported reconstructions and solvers. An explicit nonnegative
`difmag` also works with PCM/PLM; zero disables the correction. Nonfinite values and
negative values other than `-1` are rejected. The resolved coefficient is bound
to checkpoint context for PPM and positive PCM/PLM settings, so automatic and explicit equivalent values can
resume each other, while changed coefficients are rejected. PPM checkpoints
from before this context addition are incompatible. Older zero-policy PCM/PLM
checkpoints require explicit `difmag=0`. See the
[reference-default checks](validation/spray-default-difmag-20260925.md). See [PPM checks](validation/difmag-spray-input-20260925.md) and
[PLM input/restart checks](validation/spray-plm-difmag-input-20260925.md).
Regular MPI spray also accepts `difmag` in `&mpi_spray`, checks collective
agreement, and binds the resolved coefficient to restart context
([checks](validation/difmag-mpi-spray-20260925.md)). MPI AMR spray also reads the
coefficient from the gas namelist and includes
its resolved value in the shared gas checkpoint context
([checks](validation/difmag-amr-input-20260925.md)).

## Build and run

A fixed H2/O2 gas mechanism with water droplets requires no new dependencies:

```sh
cmake -S . -B build/spray -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_TESTS=ON -DPELEF_ENABLE_SPRAY=ON
cmake --build build/spray --parallel 4
ctest --test-dir build/spray -L spray --output-on-failure --no-tests=error
mkdir water-run && cd water-run
../build/spray/pelef_spray_3d ../cases/spray_3d/water.nml
```

For liquid-fuel combustion, the reference fixture is FFCM1_Red (21 species,
124 reactions) from the **already pinned** PelePhysics revision. The helper
checks the upstream YAML and license bytes before importing. It requires
Cantera 3.2.0 and never overwrites different local source data.

```sh
python3 -m pip install cantera==3.2.0
python3 tools/prepare_spray_reference.py --output-dir build/methanol-reference
cmake -S . -B build/spray-cvode -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_TESTS=ON -DPELEF_ENABLE_SPRAY=ON \
  -DPELEF_ENABLE_SUNDIALS=ON -DPELEF_ENABLE_CANTERA_REFERENCE=ON \
  -DCMAKE_PREFIX_PATH=/absolute/path/to/your/sundials/install \
  -DPELEF_MECHANISM_BUNDLE="$PWD/build/methanol-reference/ffcm1_red.json" \
  -DPELEF_MECHANISM_SOURCE="$PWD/build/methanol-reference/ffcm1_red.yaml"
cmake --build build/spray-cvode --parallel 4
ctest --test-dir build/spray-cvode -L spray --output-on-failure --no-tests=error
mkdir methanol-run && cd methanol-run
../build/spray-cvode/pelef_spray_3d_selected ../cases/spray_3d/methanol_burn.nml
```

SUNDIALS must provide the static double-precision official Fortran CVODE,
serial-vector, dense-matrix and dense-linear-solver interfaces built with a
compatible Fortran compiler. The qualification workflow pins SUNDIALS 7.2.0 at
`71a4cc9ad5e7bc8b4e33a1ca9795b4e96883f9a6`. `chemistry_method='cvode'` fails
explicitly when that backend was not built. With spray disabled (the default),
the preceding application/install target inventory remains unchanged.

Offline reference preparation uses both `--yaml PATH --license PATH` with the
same exact pinned checks. Source/license hashes and the illustrative liquid
parameters are recorded in [the manifest](../references/spray_methanol.json).
The reference helper retains the upstream license; it does not select a license
for this project.

## What the supplied cases establish

| Input or test | Purpose | Not established |
| --- | --- | --- |
| `water.nml` | Water evaporation and two-way gas feedback without reaction | Diesel fuel properties or injector accuracy |
| `methanol_burn.nml` | Closed uniform liquid-fuel evaporation through oxidation and net heat release | A spatial spray flame, lift-off, or measured ignition delay |
| `methanol_cone.nml` | Nonuniform cone, scheduled injection, transport, reacting cells and LES together | A physical inflow/nozzle boundary or turbulence statistics |
| `check_spray_cantera.py` gas-only cases | Independent full temperature history and final composition against Cantera | An external reference for the liquid model |
| `check_spray_cone.py` | All-phase conservation, elements and exact coupled restart with injection/chemistry/LES | Experimental spray validation or parallel scalability |

For a new run, output paths and checkpoint destinations must not exist. Outputs
are `.history.csv`, `.gas.csv` and `.parcels.csv`. History records gas/liquid
mass, system momentum and total energy, injection mass, balance residuals,
temperature bounds, live parcels and Sauter mean diameter. A parcel's `mass`
is the mass of **one** droplet; its multiplicity multiplies all transfer terms.

## Models and conventions

SI units are used throughout. Each liquid population is dilute, spherical, single-component,
spatially isothermal, with constant or temperature-dependent polynomial density. `vapor_name` must identify a gas
species in the selected mechanism. Liquid property numbers in the examples
are illustrative constants, **not** a fitted methanol/water property database.

The reference energy is `e_l(Tref) = h_v(Tref) - Lref`, using the selected
mechanism's NASA7 gas enthalpy. Thereafter `e_l(T) = e_l(Tref)+integral(cp_l(T), Tref, T)` and
the latent heat in the droplet energy equation is `h_v(T)-e_l(T)`. The
incompressible liquid pressure-work correction and surface energy are neglected.
This choice aligns chemical formation-energy references across phases. The
gas receives minus the parcel's discrete change of mass, momentum and **total**
energy; latent heat, evaporated kinetic energy and drag work are not added a
second time.

The surface saturation pressure uses a **constant-Lref** Clausius--Clapeyron
law anchored at `(Tref, pref)` by default. States with `psat >= pgas` are rejected
by default; the optional surface cap below limits the vapor-state evaluation temperature. Condensation is not modeled;
a saturated or supersaturated ambient produces zero evaporation. Film state
uses the one-third rule, NASA7 thermodynamics and existing molecular transport.

The vapor/carrier diffusivity is the carrier-normalized harmonic binary value,
`D_vc = sum(j!=v, X_j) / sum(j!=v, X_j/D_vj)`. It is evaluated directly from the
binary matrix and carrier composition, recovering the exact binary coefficient
for a single carrier. Do **not** apply a PelePhysics `rho*D` conversion to
PeleF's differently defined mixture-averaged coefficient. A regression reproduces
the erroneous molecular-weight conversion and guards the corrected binary limit.

Drag uses the [Schiller--Naumann sphere correlation](https://cpp.openfoam.org/dev/SchillerNaumann_8C_source.html) (`Cd*Re=24*(1+0.15*Re^0.687)`
below Re=1000; `Cd=0.44` above). Transfer supports Ranz--Marshall or an
Abramzon--Sirignano-style blowing correction. The latter uses the pinned direct
B_T iteration with absolute tolerance 1e-7 and fuel-only skin heat capacity
(supplied-film zero-field fallback uses film Cp), and fails on nonconvergence.
See the [iteration checks](validation/heat-iteration-20260925.md),
[PeleMP equations](https://amrex-combustion.github.io/PeleMP/Equations.html) and
[Abramzon and Sirignano (1989)](https://doi.org/10.1016/0017-9310(89)90043-4).
The constant-latent saturation law, drag correlation and regular-grid CIC
source distribution differ from PeleMP: this is not bitwise PeleMP parity.

Frozen-film drag and heat relaxation are integrated exponentially, evaporation
uses a nonnegative diameter-squared update, and parcel steps are limited by
relative changes, at most a small temperature change, and cell traversal.
A minimum-diameter event transfers the **entire** residual mass/momentum/energy.
By default particle locations wrap periodically. Gas interpolation and source deposition
use the same eight cloud-in-cell weights. Inadmissible gas or liquid candidates
are rejected with step reduction; exhausting the retry/substep limit rolls back
the entire caller-visible coupled step. Sequential parcel feedback is
deterministic but not parcel-order-independent or MPI distributed.

The SGS model uses `mu_t=rho*(Cs*Delta)^2*sqrt(2*S:S)`, deviatoric stress,
turbulent heat conduction and corrected species diffusion with zero net species
mass flux. Conservative face energy flux includes stress work and species
enthalpy diffusion. Its SSPRK2 update checks positivity/EOS and a diffusive
step bound. `Cs=0` is an exact no-op. There is no isotropic SGS kinetic-energy
model, wall damping, dynamic coefficient, stochastic parcel dispersion or
unresolved turbulence--chemistry interaction closure. See
[the PeleC LES discussion](https://amrex-combustion.github.io/PeleC/LES.html) for
context; this smaller closure is independently implemented and is not all of
PeleC's compressible LES model.

A coupled step applies half spray/LES/chemistry updates, the existing gas
hydro/molecular step, and reverse half updates. Each CVODE cell call uses that
cell's current density/internal energy and releases its context even on failure.
It does not reuse an integrator with stale density after evaporation. All gas,
temperature and parcel changes are committed together. Symmetric operator
ordering does **not** make the frozen-film parcel scheme second-order; temporal
and mesh convergence must be checked for each research problem.
`check_spray_convergence.py` checks history-error reduction at three outer
timesteps for the closed methanol example; it does not certify all meshes or
second-order spray accuracy.

## Injection and restart

The cone uses deterministic equal-solid-angle sampling with an azimuth based on
the monotonically assigned parcel ID. Droplets share the input temperature. `diameter_shape=0` (default) is
monodisperse at `drop_diameter`. A positive shape in [0.25,20] enables a
truncated Weibull/Rosin--Rammler **mass CDF** between `diameter_min` and
`diameter_max`, with `drop_diameter` as the scale (not the Sauter mean).
Midpoint inverse-CDF quantiles each carry equal parcel mass; multiplicity
therefore varies as diameter to the minus third power. This reproduces the
specified mass distribution, not a number-based Weibull distribution.
It is not OpenFOAM's differently interpreted `massRosinRammler` correction.
[The underlying truncated distribution](https://api.openfoam.com/2412/RosinRammler_8H_source.html)
is sampled deterministically; it is not a turbulent dispersion model. Prescribed mass flow is introduced in timed pulses;
the final pulse uses the remaining window duration so injected mass is exact.
The fluid time step is clipped to injection events. This is a prescribed
particle source in a periodic box, **not** a resolved nozzle or open inlet.

`stop_after_steps` plus `checkpoint_file` writes a coupled checkpoint;
`restart_file` resumes into a **new** output prefix. The file contains gas and
liquid state, time/step/event counters and original/injected budgets. It embeds
the mechanism identity, application version, backend availability and physical
input context; changing physics or timestep controls is refused. Output paths,
stop controls and a later final time may change. Ordered IDs, shape, finite
values, state/temperature consistency, event counts and balance ledgers are
checked before changing caller state or opening output files.

Publication closes a temporary file before an atomic no-overwrite POSIX hard
link in the destination directory. Existing checkpoints are never replaced.
This requires a Linux/POSIX filesystem supporting hard links; it is not an
`fsync`-durable guarantee against power loss. Files are strict formatted records,
not cryptographically authenticated data. Malformed/truncated/extra/null records
and changed physical context are rejected. Read failure leaves existing state
unchanged. Restart compatibility is intentionally limited to the same recorded
solver version/physics, not an indefinite checkpoint ABI.

## Explicit remaining work

Production use still needs calibrated temperature-dependent liquid properties,
boiling/real-fluid handling, multicomponent liquids, additional breakup/dispersion models and collision, calibrated size distributions
and turbulent dispersion, physical injector and wall-impact models, particle
migration with MPI/AMR/EB, scalable output/checkpointing, detailed diesel-fuel
mechanisms and performance tests. Selected chemistry uses a 2--128-species
limit and currently supported reaction forms; FFCM1_Red does not represent
HXN/HMN/AMN. No claim of measured ignition, penetration or liquid length is made.


## Regular-grid wall and outflow boundaries

`boundary_faces(6)` in `&spray` uses x-low, x-high, y-low, y-high, z-low,
z-high order. The default remains periodic. Supported alternatives are
`inflow`, `outflow`, `slip_wall`, `no_slip_wall`, and `isothermal_wall` for gas; parcels
reflect specularly at each wall type. Opposing periodic faces must be paired.
This ideal reflection is not a sticking, splash, film or impact model.

Physical faces use a clamped, normalized interpolation/deposition stencil with
unique donor cells; particles never sample the opposite side of a physical wall.
Outflow removes a parcel's remaining mass, momentum and energy from the system;
wall reflection transfers normal momentum without changing parcel energy.
The split step accumulates gas surface fluxes and parcel exchange in one ledger.
Physical-boundary history adds `boundary_mass`, `boundary_px`, `boundary_py`,
`boundary_pz`, `boundary_energy`; positive values enter the simulated system.
The reported conservation errors subtract this ledger and injected inventory.

`wall_reflection.nml` supplies an executable example. Change its faces to
`outflow` for particle escape. Checkpoint V2 preserves boundary exchange, and
its context binds all six face types. Periodic runs retain checkpoint V1 and
their existing history columns. Restart reproduces gas/parcel CSVs and the final
history row exactly in the focused [boundary checks](validation/spray-boundaries-20260915.md).

Nonperiodic faces also support the deviatoric Smagorinsky model. Physical-face
stress, heat and species fluxes use mirrored/clamped ghost states; the surface
ledger includes both LES RK stages in both split half-steps. The diffusion bound
uses only faces touching physical cells, excluding the artificial outer ghost
boundary. Set `smagorinsky_constant=.17` and `shear_velocity=10` in the example
to exercise nonzero SGS. See [focused SGS evidence](validation/sgs-boundaries-20260915.md).

This uses the existing constant-coefficient closure without wall damping or a
wall function. Near-wall turbulence qualification, characteristic inlet models,
AMR/EB/MPI parcels and physical wall-impact models
remain open. Injection is still an internal parcel source.


## ETAB representative-drop breakup

Set `breakup_model='etab'`, `surface_tension` (N/m), and `liquid_viscosity`
(Pa s); both material properties must be positive. The default is `none`.
`etab.nml` exercises constant illustrative water properties at high gas density
and relative velocity. These values are a numerical fixture, not a fitted fuel
property model.

The deformation oscillator and ETAB child-size law follow the pinned
[PelePhysics TABBreakup.H](https://github.com/Pele-Suite/PelePhysics/blob/6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b/Source/Spray/BreakupSplash/TABBreakup.H).
Breakup decreases a representative droplet's mass and increases its multiplicity
at constant parcel inventory, velocity and temperature. It preserves system
mass, momentum and the existing thermal/kinetic energy convention, which omits
surface energy. Child radii are bounded by `minimum_diameter/2`; an excessive
substep count or invalid state rejects the complete coupling step.

`deformation` and `deformation_rate` are included in breakup-enabled parcel CSVs
and checkpoint V3 (periodic) or V4 (with physical-boundary inventory). Material
parameters and model identity are bound into checkpoint context. Models without
breakup retain V1/V2. The damped oscillator advances only the remaining interval
at each internal step; zero-slip and overdamped states avoid singular formulas.

The representative-size/deformation update can also partition computational
parcels as described below. Stochastic transverse child velocities, primary atomization,
collision/coalescence and experimental breakup qualification remain open. Source hashes,
license and deliberate scope differences are recorded in
[the reference manifest](../references/spray_etab.json) and
[focused evidence](validation/etab-20260915.md).


## Mixed-liquid thermodynamic API

`spray_mixture_mod` provides the shared thermodynamic operations for up to 16
components of one liquid droplet, distinct from the existing per-port pure liquids.
Multiple components may name the same gas species. Their liquid properties and
component masses remain separate, while surface vapor and gas mass sources are
summed by gas species. Independent liquid molecular weights are supported as
described below. Ambient vapor partitioning follows the pinned `Drag.H` gate,
which compares pure-component saturation pressure with the allocated ambient
partial pressure before applying Raoult weighting. This changes humid/dilute
mixture behavior from the former gate. The updated mixture evolution identity
rejects continuation from the former model. See
[mapping evidence](validation/shared-vapor-mapping-20260925.md).

Use `liquid_molecular_weight`, `port_molecular_weights(port)` or
`mixture_molecular_weights(component)` in regular spray input, in kg/kmol.
Zero uses the mapped gas species' molecular weight; a nonzero override must
be finite and within 1e-3–1e5 kg/kmol. Zero in a port entry also means the gas
value, independently of the scalar override. MPI AMR's `amr_spray_properties`
and both EB drivers' `eb_spray` groups accept `mixture_molecular_weights`.
Unused array entries must be zero. The settings participate in restart context.

Liquid mole fractions, Clausius pressure and pressure-dependent boiling use
the liquid molecular weight. Surface vapor conversion retains the mapped gas
species' weight in its numerator, as in the reference. Gas thermodynamics,
conservative mass exchange and liquid energy references keep their existing
definitions. This is a surrogate mapping, not a declaration or validation of
independent liquid elemental formulas. The `pelemp-cap` boiling calculation
also supplies the wall-impact temperature estimate. See
[molecular-weight integration](validation/liquid-molecular-weights-20260925.md).
It accepts pure-component property records and normalized mass fractions.
Density uses additive specific volume, `1/rho = sum(Y_i/rho_i)`; heat capacity
and stored liquid energy use `sum(Y_i*cp_i)` and `sum(Y_i*e_i)`. Each component's
energy retains its mechanism formation reference. No excess mixing energy is added.

`liquid_mixture_properties` evaluates these quantities at a common temperature.
`liquid_mixture_temperature` inverts the energy with bracketed Newton iteration,
inside the active components' common density/heat-capacity fit ranges. Zero-mass
components do not constrain the temperature range. Non-normalized fractions,
duplicate component species, invalid properties or unreachable energies reject
without replacing the outputs.

`advance_mixture_heating` applies a supplied parcel heat amount at fixed composition,
recovers temperature and returns the opposite actual energy increment for the gas.
Failure preserves the parcel and returns zero transfer. The transport update
below supplies mixture heat transfer and component evaporation.
The CLI supports shared or per-port mixed-liquid compositions using the inputs below. The
mixture grid API retains component masses and supports their checkpoint persistence
and combined hydro/chemistry stepping. The
reference's optional boiling-temperature property caps are not applied here.
See [source](../references/spray_mixture.json) and
[focused evidence](validation/spray-mixture-thermo-20260916.md).

### Mixed-liquid surface and component exchange

`liquid_mixture_surface` converts mass fractions to liquid mole fractions and
returns equilibrium surface fractions `X_liquid_i * psat_i / P` (Raoult's law).
It rejects a saturated combined surface and leaves outputs unchanged on failure.
A component's pure saturation pressure may exceed ambient pressure when dilution
keeps its partial pressure below ambient. Optional per-component surface caps
use the actual total gas pressure before applying the liquid mole fraction.
This surface routine returns equilibrium fractions; the film builder below adds
ambient-vapor suppression and transport coefficients.

`advance_mixture_evaporation` takes a parcel, component masses **per physical
droplet**, requested evaporated masses on the same basis, component vapor
enthalpies in J/kg, and heat supplied by the gas to the **whole parcel** in J.
It updates the supplied component-mass array and total droplet mass, recovers the
remaining mixture temperature, and returns gas conserved-inventory increments
(mass, three momenta, total energy, followed by gas species masses).
Components map through their pure-liquid vapor indices; duplicate mappings reject.

For a surviving droplet, the target liquid energy equals old liquid energy plus
supplied heat minus evaporated vapor enthalpy. The gas receives the actual
opposite discrete liquid inventory change, including evaporated kinetic energy;
no separate latent-heat term should be applied again. Velocity and position do
not change in this exchange operation. A completely evaporated droplet returns
its finite entire inventory once, following the existing terminal-drop convention;
there is no residual liquid temperature to solve or additional heat debit.
Overdraw, invalid composition, unsupported energy states and invalid species
indices reject without modifying parcel/component masses and return zero exchange.
This API consumes calculated evaporation amounts; the frozen-gas update below
supplies them. The grid API below retains component masses alongside the parcel array.
See [component-exchange evidence](validation/spray-mixture-exchange-20260916.md).

### Mixed-liquid film transport and frozen-gas update

`build_mixture_film` combines liquid surface fractions with the ambient gas.
The pinned pure-saturation-pressure gate described above selects vapor components.
The complete ambient mixture, including existing fuel vapor, supplies the mean
molecular weight in the vapor-to-mass conversion. The one-third film rule supplies
droplet skin temperature and composition; wall film uses its one-half rule.
Inactive gas species are renormalized to close the skin mass fractions. Single
liquid and mixture builders use the same conversion and activity ordering.
Each vapor uses the
existing carrier-mole harmonic binary diffusivity. Effective diffusivity is the
surface-mole-weighted average; component evaporation shares are proportional to
`X_surface_i * D_i`, following the reference's component weighting structure.

`mixture_transfer_rates` applies the existing selectable drag and
Ranz–Marshall/Abramzon–Sirignano correlations using the actual mixed-liquid density
and heat capacity. It returns inverse drag time, heat coefficient in W/K,
component mass rates in kg/s **per physical droplet**, and vapor enthalpies.
For Abramzon–Sirignano heat correction, the builders supply skin-temperature Cp
averaged over mapped fuel gas species only, including inactive mapped species
present in the skin. Aliases of one gas species are counted once. Manually supplied
`spray_film` values retain the skin-Cp fallback when `vapor_cp=0`.
The mass Spalding number is capped at 20 for Abramzon–Sirignano and wall-film
evaporation; the selectable Ranz–Marshall path is unchanged. PeleF retains its
binary diffusivity convention, guarded small-B limit and nonconvergence rejection;
the direct B_T iteration follows the reference for ordinary positive B_M.
This is not full upstream coefficient or integrator parity. See
[fuel-Cp evidence](validation/vapor-cp-20260925.md).
The ambient-mixture conversion
and wet-gas restart checks are recorded in
[the skin-state evidence](validation/ambient-vapor-skin-20260925.md). Regular
checkpoint context and mixture evolution identity changed with this numerical
correction, so older model contexts are rejected.

`advance_mixture_frozen_gas` combines exponential drag, trapezoidal flight,
heat transfer and component evaporation over a supplied interval. Gas state is
held fixed during this operation; returned conserved increments include drag
work and species exchange. It changes the parcel and external component-mass
array only on success. Overdraw or inadmissible temperatures reject the interval;
the caller must shorten it. Evaporating trace components at or below a mass
fraction of 1e-12, or with a predicted residual below that fraction of the current
droplet mass, release their entire residual inventory. They no longer constrain
the grid step through a relative component-mass limit. Suppressed components are
retained. If evaporation reduces the predicted mixture diameter below
`minimum_diameter`, every component is transferred to the gas using the terminal
finite-inventory convention. These transfers conserve mass and energy rather
than discarding small species and renormalizing their survivors.
A zero interval preserves state exactly. This remains
a composable API used by the grid and wall coupling below.
See [rate and update evidence](validation/spray-mixture-rates-20260916.md).

### Mixed-liquid grid coupling and composition restart

`advance_mixture_coupling_3d` advances a parcel array and a persistent
`component_mass(component,parcel)` array together with a regular 3D gas grid.
All parcels use the supplied ordered component table. Each column stores masses
per physical droplet and must sum to that parcel's droplet mass. The routine
samples gas with the existing particle stencil, limits drag/temperature/component
changes and particle travel, and halves rejected feedback intervals. Gas,
particle, component and boundary states commit only after the full update succeeds.

The boundary API supports periodic wrapping, specular wall reflection and
outflow/inflow removal. Boundary budgets include each component's species mass;
`mixture_system_totals` includes both gas and liquid inventories on the same basis.
ETAB/KHRT representative-size updates preserve composition while scaling droplet
mass and multiplicity. Their surface tension and liquid viscosity must be supplied
as mixture properties through the existing options. Computational splitting copies
component masses per droplet into its descendants. Optional KH shed-mass children
inherit the parent composition, with conservative transverse kinetic energy debited
from gas. Children advance through the remaining step; capacity or gas-state failure
rolls back the entire update. KH accumulators and composition survive computational
splitting and restart. Wall-impact/film conversion is available through the optional coupled arguments below.
See [mixed KH evidence](validation/spray-mixture-kh-20260916.md).

Both checkpoint routines accept optional `component_mass`. The `_MIX_V1` suffix
adds a component count after grid dimensions and a component-mass row after each
parcel record. Counts, shapes, nonnegative finite masses and parcel-mass closure
are checked. A reader without composition output rejects this format; invalid
reads preserve existing state. The caller's restart context must include `mixture_evolution_id` and identify the
ordered component properties and mixture settings, as with existing pure-liquid
checkpoint contexts. The reader commits the allocated component matrix together
with the gas and parcel state.

These regular-grid coupling/checkpoint APIs are used by the mixture namelist
and scheduled injection path below.
See [grid/restart evidence](validation/spray-mixture-grid-20260916.md).

### Mixed-liquid reacting-flow step

`advance_spray_reactive_3d` accepts the optional pair `mixture_components` and
`component_mass`. It advances the same mixture state in both spray half-steps,
with the existing hydro, molecular transport, LES and chemistry operators between
them. External acceleration retains its impulse/work ledger. Gas, particles and
component masses commit together only after every stage succeeds; late chemistry
or transport rejection returns zero external/boundary ledgers and preserves inputs.

The mixture pair cannot be combined with the per-port pure-liquid property tables.
Mixed particle wall impact requires `wall_films` and `film_component_mass` together.
Gas wall conditions remain available. No new chemistry
backend or time-order accuracy claim is introduced. The focused check exercises
the native implicit backend, not a new CVODE qualification. See
[reacting-mixture evidence](validation/spray-mixture-reactive-20260916.md).

### Mixed-liquid application inputs

Set `mixture_count` to 1–16 and supply that many entries in:

| Input | Meaning / units |
| --- | --- |
| `mixture_vapor_names` | Distinct gas species names |
| `mixture_fractions` | Initial liquid mass fractions, nonnegative and summing to one |
| `mixture_densities` | Constant-model component densities, kg/m³ |
| `mixture_cps` | Constant-model component heat capacities, J/(kg K) |
| `mixture_reference_temperatures` | Component reference temperatures, K |
| `mixture_reference_pressures` | Reference saturation pressures, Pa; default 101325 |
| `mixture_latent_heats` | Reference latent heats, J/kg |

These inputs define the component formation-referenced energy and Clausius vapor
pressure. Unused entries must remain at their defaults. `mixture_count=0` keeps
the existing pure-liquid mode. Component property models can be set independently:

| Input | Values / layout |
| --- | --- |
| `mixture_density_models`, `mixture_cp_models` | `constant` (default) or `polynomial` |
| `mixture_density_coefficients`, `mixture_cp_coefficients` | `(4,16)` cubic coefficients, `a + T*(b + T*(c + T*d))`, SI units, T in K |
| `mixture_density_temperature_ranges`, `mixture_cp_temperature_ranges` | `(2,16)` valid temperature bounds in K |
| `mixture_vapor_pressure_models` | `clausius` (default) or `antoine` |
| `mixture_antoine_coefficients` | `(4,16)` values A,B,C,D for `psat = D*10**(A-B/(T+C))`, Pa and K |
| `mixture_antoine_temperature_ranges` | `(2,16)` valid temperature bounds in K |
| `mixture_boiling_models` | `reject` (default) or `surface-cap`, the latter requires Clausius |
| `mixture_critical_temperatures` | Component critical temperatures in K for surface-cap |

For example, `mixture_cp_coefficients(:,1)=4000,1,0,0` selects the first
component's coefficient column. Coefficients and ranges for inactive models must
remain zero, and unused component entries must retain defaults. Component energy
uses the integrated heat-capacity polynomial relative to its reference temperature.
The surface cap limits the pure-component vapor-state temperature before Raoult
partitioning; it does not implement physical mixture boiling. Models, coefficients
and bounds enter restart identity. Default constant-property contexts are preserved.
See [component property evidence](validation/spray-mixture-properties-20260916.md).

Initial and scheduled parcels use `mixture_fractions` by default. To override a
port, set `port_mixture_fractions(component,port)`, a `(16,16)` array. For example,
`port_mixture_fractions(1:2,1)=.2,.8` assigns the first port a 20/80 mass blend.
Each nonzero column must be nonnegative and sum to one; an all-zero column inherits
the shared composition. Rows beyond `mixture_count` and columns beyond
`injection_port_count` must stay zero. Components may be absent from individual
ports; every port uses the same ordered component property table. Keep a valid
shared composition even when overriding all active ports.

Both initial and scheduled injections use each port's selected composition for
density, energy, component inventories and rate-derived velocity. Restart identity
includes active port compositions. See [port-composition evidence](validation/spray-port-mixture-20260916.md).
Existing size laws,
nozzle geometry, per-port temperatures and shared/independent schedules remain
available. Injection density and energy use the actual mixture at each port's
injection temperature; rate-derived velocity uses that density. New component
masses append alongside each injected parcel, with species-resolved mass ledgers.
Per-port pure-liquid property overrides cannot be combined with mixture mode. Prescribed mixture viscosity/surface tension permit representative
ETAB/KHRT, KH shed children and computational splitting.

Parcel CSVs add `mass_<species>` columns in kg per physical droplet. Their sum is
`mass`; `diameter` and Sauter mean diameter use the evolving mixture density.
The ordered component table and initial fractions enter restart compatibility.
Changed composition/properties reject restart, and checkpoints retain every
parcel's evolved component masses. The terminal/trace update changes the
mixture evolution tag to `mixed-liquid-components-v2-terminal`; earlier mixture
checkpoints reject rather than silently changing their continuation policy.
See [terminal-state evidence](validation/spray-mixture-terminal-20260916.md).

[mixture.nml](../cases/spray_3d/mixture.nml) is a small synthetic example. The
illustrative H2O/H2O2 properties demonstrate software behavior, not a validated
liquid blend. See [CLI/injection evidence](validation/spray-mixture-cli-20260916.md).

### Mixed-liquid wall-impact event API

`resolve_mixture_wall_impact` partitions an incident mixture among the retained
parent, splash children and a deposited film patch using the existing six-regime
impact kernel. Parent and child component masses are kg per physical drop;
`film_mass` is kg per complete film patch. All outgoing material inherits the
incident composition and temperature. The kernel evaluates mixture density and
energy at that temperature and returns wall momentum/energy exchange. Failed
classification or child/film construction preserves all supplied state arrays.

The caller supplies the dimensionless impact parameters, including wall-to-boiling
temperature ratio. The planar dispatcher below evaluates that ratio. Grid/CLI coupling remains to be connected. The checkpoint API below retains
mixture-film composition. The CLI supports this mixed-liquid wall path. There is no entrainment
from an existing film. See [event evidence](validation/spray-mixture-wall-event-20260916.md).

### Mixed-liquid film exchange API

`advance_mixture_wall_film_exchange` accepts component evaporation amounts in kg,
component vapor enthalpies in J/kg and gas/wall heat amounts in J into a stationary
film patch. It inverts the remaining mixture energy for temperature and updates
height at fixed footprint using the new density. `component_mass` stores total kg
for the patch. Returned gas changes include individual vapor species; `wall_energy`
is the opposite of supplied wall heat. Momentum exchange is zero for this stationary
film operation.

Complete evaporation and an evaporation-driven height below `minimum_height` return
all remaining species and energy to gas. The full inventory convention also applies
when the final temperature cannot be defined because no liquid remains. Invalid
component loss or temperature inversion leaves the patch and composition untouched.
This supplied-exchange API is used by the frozen-rate update below. Mixed films
are not yet connected to the grid driver. See
[film exchange evidence](validation/spray-mixture-film-exchange-20260916.md).

### Mixed-liquid film rates and frozen update

`mixture_wall_film_rates` evaluates Raoult surface fractions and diffusion-weighted
component evaporation through the existing stationary-film correlation. Gas-side
properties use the wall half-composition rule and
`Tskin = Tfilm + 0.5*max(0,Tgas-Tfilm)`. Gas and wall heat rates use the film thickness
and supplied liquid conductivity. The temperature-rate estimate subtracts each
component's evaporation rate times its own vapor-minus-liquid enthalpy, divided by
the film mass and mixture heat capacity.

`advance_mixture_wall_film_frozen` holds these rates fixed over a supplied interval
and calls the conservative exchange API. Actively evaporating components whose
current or predicted residual mass is at most `1e-12` of film mass are returned
entirely to gas. Ambient-suppressed components stay in the film; a zero interval
changes nothing. The grid depletion limiter excludes already trace components.
Bulk component overdraw or an
inadmissible remaining temperature rejects; a grid driver must reduce the step in those
cases. It retains the full-inventory minimum-height removal policy. Composition
transport, thermal updates and rates are now available as APIs; grid/CLI feedback
and application restart wiring remain open. See
[rate evidence](validation/spray-mixture-film-rates-20260916.md).

### Mixed-liquid planar wall dispatch

`mixture_impact_at_planar_wall` finds the earliest crossed box face, calculates
impact parameters with the current mixture density and prescribed mixture viscosity
and surface tension, and dispatches the existing impact event. Outgoing parcel
component arrays use kg per drop; deposition components use whole-patch kg.
No crossing leaves the supplied state unchanged; failed dispatch rolls back.

`mixture_impact_boiling_temperature` uses the mass-weighted pure-component boiling
temperatures for the wall-temperature ratio, following pinned `boilT_mix` and
AhamedSplash. Component temperatures use the existing Clausius or bounded Antoine
inversion at local pressure. This empirical impact temperature is not a Raoult
bubble-point calculation or a physical boiling-rate model. Existing film geometry
supplies local thickness; film entrainment remains absent. See
[crossing evidence](validation/spray-mixture-wall-crossing-20260916.md).

### Mixed-film checkpoint API

The checkpoint reader/writer accept `film_component_mass(component,patch)` in kg
per complete film patch. It accompanies `component_mass`, `film_state=.true.` and
`wall_films`; the format adds `_FMIX_V1` and a component row after each film record.
The number and order of components match the particle component matrix. Callers
must bind the ordered component property table and film settings in the context.

A mixed checkpoint with film state requires this matrix, including a zero-column
matrix when there are no films. Shape mismatch, negative/nonfinite component mass
or failure to sum to film mass rejects. Reading validates the complete candidate
before changing gas, particles, films or either component matrix. Existing pure-film
and particle-only mixed formats retain their layout. The CLI retains this mixed-film checkpoint state.
See [checkpoint evidence](validation/spray-mixture-film-checkpoint-20260916.md).

### Mixed-film grid feedback API

`advance_mixture_film_coupling_3d` advances stationary film patches and their
component matrix against a regular gas grid. It samples gas properties at each
patch, limits substeps by component depletion and temperature change, and deposits
component mass and energy with the particle interpolation weights. Failed gas
recovery halves the attempted interval; failure of the complete update preserves
gas, temperatures, films and component masses. Optional `start_times` delays each
patch until its birth time within the interval.

The returned boundary ledger records heat entering the simulated gas/film system
from walls. Fully evaporated patches and corresponding matrix columns are removed
together. `mixture_system_totals` accepts optional `wall_films` and `film_mass` to
include stationary-film mass, species and internal energy. This API is connected to particle impact and reactive/CLI stepping. See [grid evidence](validation/spray-mixture-film-grid-20260916.md).

### Combined mixed-particle and wall-film update

`advance_mixture_coupling_3d` accepts optional `wall_options`, `wall_films` and
`film_component_mass`. The two film arrays must accompany each other. Enabled wall
impact requires these arrays, even when initially empty. It processes the accepted
particle flight's first wall crossing before aerodynamic breakup; a hit skips
breakup for that substep. Splash children and deposited patches inherit the incident
composition and start at the accepted substep's end.

After particle flight, stationary films exchange mass and energy with the grid for
their remaining interval. Particle and film substeps share the configured budget.
Capacity failure, inadmissible impact or rejected film feedback rolls back the entire
particle/gas/film update. Boundary ledgers combine wall impulse/energy and wall heat.
Reactive-driver and CLI integration are available below. See
[combined evidence](validation/spray-mixture-wall-coupled-20260916.md).

### Reactive mixed-wall integration

`advance_spray_reactive_3d` accepts `film_component_mass` alongside mixed-particle
composition and `wall_films`. Both spray half-steps now carry wall impact and mixed
film evolution through the same candidate state as hydro, molecular transport,
LES and chemistry. The complete state commits only after every stage succeeds;
a late chemistry failure restores both component matrices with gas, particles and
films. Pure-material property tables still cannot accompany mixture tables.

The existing prohibition on moving particle-impact walls remains. CLI input/output
and checkpoint wiring retain both particle and film compositions. See
[reactive evidence](validation/spray-mixture-wall-reactive-20260916.md).

Mixed-film evolution now uses `mixture_film_evolution_id='mixed-film-v2-trace'`
in CLI checkpoint context. Earlier mixed-wall checkpoints reject because their
trace-depletion evolution differs. The storage layout, pure-wall and particle-only
mixture contexts are unchanged. API callers must include the evolution ID in their
own context. See [trace evidence](validation/spray-mixture-film-trace-20260916.md).

### Mixed-wall CLI and outputs

Combine mixture inputs with the existing `wall_impact_faces`, `surface_tension`,
`liquid_viscosity`, `wall_film_evolution` and `wall_liquid_conductivity` settings.
Stationary film heat conductivity is prescribed; this is not a calibrated mixture
conductivity model. The `.films.csv` output adds `mass_<species>` columns in kg per
whole patch, while parcel component columns remain kg per physical drop. History
and conservation ledgers include gas, particles and mixed films.

Restart retains both component matrices and validates film volume using the current
film composition and temperature. [mixture_wall.nml](../cases/spray_3d/mixture_wall.nml)
is a small synthetic deposition/evaporation example. See
[CLI evidence](validation/spray-mixture-wall-cli-20260916.md).

## EB and exterior drift API

`drift_spray_level_set` accepts the optional group `outer_faces(6)`,
`specific_energy(:)` and `outer_change(5)`. The energy array supplies each
parcel's specific liquid internal energy (J/kg). The routine processes the first
EB or root-domain crossing, then continues the remaining flight. Exterior
`outflow`/`inflow` removes represented inventory; supported wall labels reflect
specularly; paired periodic faces wrap positions. Periodic nodal geometry must
match exactly at the paired faces. The shared hit budget covers both event types.

`change` reports EB impulse/energy only; `outer_change` reports exterior system
change, including escaping liquid internal and kinetic energy. Escaped records
retain their IDs with zero mass; callers must retire their component inventories
and account for species loss. All parcels and both ledgers publish together.
Omitting the optional group retains root-contained behavior. Impact/film flight
and application exposure remain separate work. See
[focused evidence](validation/spray-eb-outer-drift-20260918.md).

The allocated serial source accepts `outer_faces(6)` with `outer_change(nvar)`
and invokes this drift after gas exchange. It computes each live parcel's liquid
energy from its component state, clears escaped component masses, and accumulates
their mapped vapor-species losses with exterior mass/momentum/energy. Its existing
five-entry boundary output continues to report EB wall exchange. Both ledgers,
gas, particles and IDs commit only after the complete interval succeeds. KH
eligibility checks stop at the first exterior crossing before querying the EB.
With impact-film arguments, all six exterior faces must be `outflow`. The
impact interval compacts escaped records and maps component losses to their vapor
species, summing components that share a species. The allocated reactive wrapper
also accepts `particle_outflow` with films through both source halves, preserving
the full-step transaction. The serial application accepts this combination and
preserves sparse surviving IDs and films through restart; see
[application evidence](validation/spray-eb-splash-outflow-application-20260918.md) and
[combined coupling evidence](validation/spray-eb-impact-outflow-coupling-20260918.md).
The allocated reactive wrapper and serial application expose `particle_outflow` for
all exterior faces, matching the existing gas outflow convention. Source losses
join the reactive boundary ledger and persist through checkpoint/restart. Other
exterior particle modes remain at the source API level. See
[application evidence](validation/spray-eb-outflow-application-20260918.md) and
[source evidence](validation/spray-eb-outer-source-20260918.md).

## KH-RT representative-drop breakup

Set `breakup_model='khrt'` with positive `surface_tension` and
`liquid_viscosity`. `khrt.nml` uses the pinned PelePhysics defaults
`khrt_b0=.61`, `khrt_b1=7`, `khrt_c3=1`. KH waves relax the representative
radius; RT instability uses a persistent waiting time, initialized to -1
and enabled by the first KH event. The Reynolds number uses film density and
viscosity; the wave equations use bulk gas density and relative velocity.
Mass reduction increases multiplicity at fixed represented mass, velocity and
temperature. The minimum diameter and transactional failure rules still apply.

The default implements the upstream `make_new_drops=false` branch. Optional
KH shed-mass children are described below. Computational splitting partitions
identical representative droplets. Surface energy is omitted.
Zero RT acceleration resets an active timer to avoid singular expressions.
Exact upstream timestep parity and physical validation have not been established.

Parcel CSVs include `rt_time`. KH-RT checkpoints use V7/V8 (periodic/physical)
or V9/V10 with computational splitting, retaining the RT timer and binding all
three coefficients into the restart context. Existing V1--V6 formats are unchanged.
See [source provenance](../references/spray_khrt.json) and
[focused evidence](validation/khrt-20260915.md).


### KH shed-mass children

Set `khrt_children=.true.` and positive `khrt_reference_diameter_cubed`
with `breakup_model='khrt'`. The reference quantity has units of m³ and is
`mean_over_jets(N_injection * mean_diameter**3)`, matching the upstream
`avg_inject_mass` normalization despite that variable's name. Supply it explicitly;
it is neither liquid mass in kg nor an automatically inferred injection dose.

The model accumulates `N_initial * (r_old**3-r_new**3)` during KH events.
This inventory remains included in the parent until release. The reference's
3% jet threshold and 90% parent threshold select partial or full release;
partial release also requires a resolved child radius and sufficient child count.
The volume denominator's numerical epsilon is converted from cm³ to m³.
The configured `minimum_diameter` supplies the radius floor, rather than the
reference's density-dependent minimum-mass rule.

Each release creates four equal-weight children in opposite transverse pairs,
with speed `0.188 * omega_KH * lambda_KH`. They inherit the parent's position,
temperature and liquid material, retain `parent_id`, and start fresh breakup states.
The transverse plane uses the parent velocity (relative gas velocity if stationary).
Gas energy decreases by the added transverse kinetic energy, closing total energy;
the paired directions preserve mean momentum. This deterministic quadrature is a
PeleF choice and does not reproduce the upstream random azimuth/count algorithm.
Creation occurs at the accepted parcel substep's end; children advance through
any remaining interval. Failed feedback, capacity checks or later chemistry roll
back the entire coupled update.

CSV fields `khrt_initial_count` and `khrt_shed_volume` expose the accumulator;
checkpoint suffix `_SHED_V1` retains it alongside RT time and parent lineage.
The enable flag and reference quantity participate in restart compatibility.
Stationary wall impact and computational splitting are supported. For splitting,
both KH accumulators are partitioned with the represented drop count, so a split
does not duplicate pending shed mass. The final piece carries rounding remainders.
The jet release threshold still applies per parcel against the prescribed global
reference; changing computational multiplicity can therefore change release timing.
When a substep crosses an enabled wall, its incident parcel undergoes wall impact
before any KH release. A hit skips aerodynamic breakup for that substep; otherwise
KH children are born at its end and can impact during their remaining flight.
This prevents the same incident mass from being released and deposited twice.
Splash children start with empty KH inventories; deposition consumes the parent's
inventory, thermal breakup resets it, and a rebound retains it. Existing wall
correlation rejection rules still apply, including unsupported wet-splash states.
See [combined wall/KH evidence](validation/khrt-wall-20260916.md).
RT still changes the representative size; it does not itself create children.
Collision, general turbulent dispersion, primary atomization and physical
qualification remain open. See [focused evidence](validation/khrt-children-20260916.md).

### Computational parcel splitting

For ETAB or KH-RT runs, set `maximum_multiplicity` to a positive upper bound on represented
droplets per computational parcel. Zero (the default) disables this partition.
`maximum_parcels` bounds allocated particle storage, up to 1,000,000.
After a coupling interval, overweight live parcels are partitioned into equal
weights, with the final child carrying the rounding remainder. Droplet mass,
position, velocity, temperature, deformation and RT timer are copied. KH initial
count and pending shed volume are divided among the pieces; aggregate mass,
momentum and energy are conserved. This increases computational parcel count
without introducing a transverse velocity or dispersion model.

Existing IDs remain fixed. New IDs append in parent order; each new child stores
the ID of the parcel it was split from. Injection roots have `parent_id=0` and
later injections continue after the last allocated ID. Restart validates parent
links and counts roots against injection events, independently of child count.
The complete split plan is checked against storage and ID limits before arrays
are replaced; failure rolls back the coupled gas/particle update.

ETAB splitting-enabled checkpoints use V5/V6, storing exact integer parent IDs in
addition to the ETAB and optional boundary state. The multiplicity/storage
limits are included in the context; parcel CSVs gain `parent_id`. For example,
add `maximum_multiplicity=1e7, maximum_parcels=1000` to `etab.nml`.
See [splitting/restart evidence](validation/spray-splitting-20260915.md) and
[KH shedding partition evidence](validation/khrt-splitting-20260916.md).


## Gas reservoir inlet

Use `boundary_faces='inflow','outflow','periodic','periodic','periodic','periodic'`
with `inflow_temperature(1)`, `inflow_pressure(1)` and `inflow_velocity(:,1)`
for an x-low gas reservoir. The six-face ordering is unchanged. The defaults are
1200 K, 101325 Pa and zero velocity; active inlet values must be finite and
thermodynamically valid. The state is prescribed in ghost cells for hydro,
molecular transport and SGS. This is not a characteristic nonreflecting inlet
or a specified numerical mass-flux boundary.

Inlet composition defaults to the initial gas mixture. Override it with
`inflow_composition_count`, `inflow_composition_species(32)` and
`inflow_composition_mole_fractions(32)`, using named species of the loaded
mechanism. A shared inlet mixture applies to all inlet faces; thermodynamic
state and velocity may differ by face. `inflow.nml` exercises an H2/O2/N2
reservoir feeding initially pure N2 with evaporating water parcels, molecular
transport and Cs=0.17.

The timestep is capped by the inlet wave speed and, when enabled, its molecular
diffusion limit. Gas convective/diffusive/SGS exchange joins the surface ledger;
all active inlet settings and resolved composition bind checkpoint context.
Particles crossing out through an inlet leave the system and are accounted as
particle outflow. A gas inlet does not create parcels; parcel injection remains
controlled by the existing separate source schedule.

See [inlet restart and conservation evidence](validation/spray-inflow-20260915.md).

Regular fixed/selected spray inputs also accept `inflow_mach(6)` (default -1).
A value from 0 to 100 enables inward normal Mach velocity on that `inflow` face;
leave `inflow_velocity(:,face)=0`. Temperature remains prescribed.
`inflow_pressure(face)=0` extrapolates current interior pressure; positive
pressure prescribes it (maximum 1e15 Pa). Mach faces use current interior
composition, including gas changed by particle feedback, rather than the shared
reservoir mixture. Hydro, molecular transport and static/Dynamic LES reconstruct
these states at each stage. Current padded inlet states also constrain each
step's wave-speed and molecular diffusion limits; the existing retry loop
handles rejected coupled candidates. For an extrapolated x-low inlet:

```fortran
inflow_mach(1)=0.2
inflow_temperature(1)=1100
inflow_pressure(1)=0
inflow_velocity(:,1)=0,0,0
```

Mach settings bind checkpoint context; changing Mach, pressure or temperature
rejects continuation. Ordinary inputs retain their previous context encoding.
The regular coupled API appends optional `inflow_mach` and `inflow_pressure`;
pass zero pressure for inactive Mach faces at that API level. The namelist
uses the existing reservoir-pressure entries and masks inactive Mach faces
before calling it. MPI/AMR spray Mach dispatch remains open.
See [Mach inlet restart evidence](validation/mach-inlet-spray-20260925.md).



## Isothermal and tangentially moving walls

`wall_temperature(6)` specifies kelvin for `isothermal_wall` faces.
`wall_velocity(3,6)` specifies a tangential velocity on `no_slip_wall` or
`isothermal_wall` faces; its normal component must be zero. Defaults are
300 K and zero velocity. Heated or moving walls require molecular transport.
`thermal_moving_wall.nml` supplies a small executable case.

Molecular wall heat transfer and viscous work use the existing regular-grid
transport operator. SGS uses the same reflected wall velocity and applies a
wall-to-cell thermal gradient over half a cell. Thermal gradients enter the
SGS timestep limit; the molecular wall limit also caps the application step.
All wall fluxes enter the accumulated momentum/energy ledger. Active wall
parameters bind restart context, so changed wall motion or temperature rejects
an incompatible restart.

Parcel impact remains specular: normal velocity reverses, tangential velocity
and temperature remain unchanged. Wall heat exchange with droplets, splash,
adhesion, roughness and wall-function calibration remain open. The gas wall
motion is tangential motion at a fixed boundary, not moving geometry.
See [focused checks](validation/spray-thermal-wall-20260915.md).


## Uniform external acceleration

`acceleration=0,0,-9.81` applies a prescribed constant acceleration in m/s² to
both gas and droplets; the default is zero. See `acceleration.nml`. Two exact
half-step velocity kicks surround the coupled operators. Gas total energy
receives the corresponding kinetic work; droplet temperature is unchanged.
The midpoint particle drift is exact for constant acceleration without drag.
The overall evaporating, coupled method retains its existing accuracy limits.

External impulse and work are separate from gas/liquid feedback and physical
boundary exchange. History gains `source_px`, `source_py`, `source_pz` and
`source_energy`; conservation residuals subtract these accumulated sources.
Checkpoints with acceleration append `_FORCE_V1` to their existing format tag
and store a source ledger. The acceleration vector binds restart context.
Zero-acceleration files retain their previous format. Step rejection rolls back
both phases and external-source increments. The timestep includes an
acceleration displacement limit, and kicks check thermal-state recovery.

This is a uniform body force, not a hydrostatic well-balanced discretization,
spatial gravity field, buoyancy closure or particle-only external force model.
Potential energy is not included in reported energy; external work accounts for
its effect. See [focused evidence](validation/spray-acceleration-20260915.md).


## Larger selected mechanisms

The importer, generator, selected input arrays, reactive state/transport checks
and MPI species-comparison storage now share a 128-species capacity. Runtime
Fortran limits use `constants_mod::mechanism_species_capacity`. Reaction and
thermo families remain restricted to the documented supported forms.
The capacity is a storage/input bound, not a performance qualification.

A 53-species, 325-reaction GRI-Mech case now runs through the existing selected
CVODE 0D and regular 3D spray executables. The 0D trajectory agrees with an
independent Cantera reactor, and the coupled chemistry/transport/spray case
restarts exactly. See [commands and evidence](validation/species-capacity-20260915.md).
MPI and AMR with more than 32 species have not been newly qualified.


## SRI pressure-dependent reactions

Selected mechanism ingestion accepts Cantera `SRI` falloff data with required
`A`, `B`, `C` and optional `D=1`, `E=0`. The normalized bundle uses
`type="sri-falloff"` and an `sri` coefficient mapping. SRI is distinct from
Troe and cannot specify both models. A and C must be nonnegative, D positive,
and A and C cannot both vanish. All coefficients must be finite.

Fortran evaluates the five-parameter broadening factor using log-sum-exp and
provides the analytical concentration derivative used by the existing kinetic
Jacobian. AMR/EB checkpoint mechanism records include the five SRI coefficients
for SRI reactions; MPI mechanism comparisons include them too. Existing
non-SRI checkpoint records keep their layout. This extends reaction kinetics,
not the liquid property model. See [focused evidence](validation/sri-20260915.md).


## Chemically activated reactions

Cantera `type: chemically-activated` reactions now use the same supported
Lindemann, Troe, SRI or Tsang broadening data as falloff reactions. Their rate is
`k0 F / (1 + Pr)`, with `Pr = k0 [M] / kinfinity`; the runtime supplies the
corresponding analytic concentration derivative. This is distinct from the
ordinary falloff prefactor. Normalized bundles use `chemically-activated`
or `sri-chemically-activated` as their reaction type.

Reaction kinds identify the different pressure laws in MPI mechanism comparison
and checkpoint records. SRI coefficients persist for both pressure laws, and a
changed coefficient rejects AMR restart. The existing reduced-pressure floor
also applies to the broadening calculation; below it the broadening derivative
is zero for the activated branch. See [focused evidence](validation/activated-rates-20260915.md).


## Tsang broadening

The selected importer accepts Cantera `Tsang: {A: ..., B: ...}` data for both
falloff and chemically activated reactions. B defaults to zero. Normalized
bundles use `tsang-falloff` or `tsang-chemically-activated` and a `tsang`
coefficient mapping. The center factor is `A+B*T`; the rest of the broadening
and its concentration derivative use the existing Troe expressions. Both
coefficients must be finite and the center must be positive at the evaluated
temperature. A nonpositive center rejects the rate calculation.

Tsang has separate reaction kinds and coefficient storage. AMR/EB checkpoint
records and MPI mechanism comparisons include both coefficients. Legacy
non-Tsang checkpoint layouts are preserved. See
[focused evidence](validation/tsang-20260915.md).


## PLOG pressure interpolation

The importer accepts Cantera `pressure-dependent-Arrhenius` reactions and emits
normalized `type="plog"` with ordered `plog` entries containing P (Pa), A, b and
Ea (J/kmol). Up to 64 entries per reaction are stored. Equal-pressure entries
remain separate so their temperature-dependent rates can be summed at runtime.
Signed A coefficients are supported when each pressure group's total is positive.
Nonpositive group totals or unrepresentable rates reject the calculation.

The runtime sums each group using scaled signed exponentials, then interpolates
log(rate) against log(pressure). Outside the supplied pressure interval it uses
the nearest endpoint rate. Pressure is R*T*sum(concentrations); its analytic
concentration derivative feeds the existing Jacobian. Internal knot derivatives
use the right interval, and endpoint branches have zero pressure derivative;
these derivatives need not be continuous at knots.

PLOG count and coefficients participate in MPI mechanism comparison and
AMR/EB checkpoint records. Existing non-PLOG checkpoint records retain their
layout. See [focused evidence](validation/plog-20260915.md).

## Chebyshev temperature/pressure rates

Selected bundles accept Chebyshev coefficient matrices up to 8 by 8, with
positive, ordered temperature and pressure fit bounds. Coefficients come from
Cantera's normalized rate API in m/kmol/s units. The runtime evaluates the
polynomial in inverse temperature and log pressure, with pressure R*T*sum(C),
and supplies its analytic derivative with respect to every concentration.
It evaluates the polynomial outside the fit interval when representable;
this extrapolation does not extend the physical validity of a mechanism fit.
Zero pressure, malformed tables and unrepresentable logarithmic rates fail.

The matrix, shape and bounds participate in MPI mechanism comparisons and
AMR/EB checkpoint contexts. Existing reaction kinds retain their checkpoint
layouts. See [the focused Chebyshev evidence](validation/chebyshev-20260915.md).

## Lu n-dodecane gas mechanism

`prepare_spray_reference.py --fuel dodecane` prepares the pinned 53-species,
268-reaction Lu mechanism using the same source/license checks as the default
methanol reference. Configure the selected bundle with its generated JSON and
original YAML. The importer and runtime preserve finite Troe A values outside
(0,1), including A=1, and signed T2; T1/T3 remain positive and the evaluated
center must be positive and representable.

The selected CVODE path now evaluates negative Newton trial concentrations
through a nonnegative extension and matching Jacobian columns. Accepted states
retain the previous composition, conservation and CVODE constraints. This
allows trace species to approach zero during stiff ignition without rejecting
the internal Newton iteration prematurely.

A phi=0.5, 1100 K, 4 MPa constant-volume ignition case reaches 2478 K at 0.6 ms
and agrees with Cantera. A separate short coupled case exercises this gas
mechanism with water parcels and exact restart; it does not supply a liquid
n-dodecane property model. See [commands and evidence](validation/dodecane-20260915.md).

## Antoine saturation-pressure fits

Set `vapor_pressure_model='antoine'`, `antoine_coefficients=A,B,C,D`, and
`antoine_temperature_range=Tmin,Tmax` to use
`psat = D * 10**(A - B/(T+C))`. T, B and C use kelvin; D scales the resulting
pressure to Pa. For coefficients tabulated in bar, D=100000. The convention
matches the pinned PelePhysics `SprayPropertiesPeleMP.H` implementation.

The default `vapor_pressure_model='clausius'` retains the existing relation
and requires the Antoine arrays to be zero. Antoine inputs require finite
coefficients, positive B/D, ordered positive bounds and T+C>0 throughout the
interval. Evaluation outside the fit range or at psat>=ambient pressure fails
transactionally. This does not implement boiling, supercritical evaporation,
multicomponent liquid properties. Temperature-dependent heat capacity is configured separately below.

`cases/spray_3d/water_antoine.nml` supplies the NIST water fit for 273--303 K
and starts droplets at 290 K. The coefficient source and units are recorded in
[`references/spray_antoine.json`](../references/spray_antoine.json); density/cp
remain illustrative constants. Changing either coefficients or bounds rejects
a checkpoint restart. See [focused evidence](validation/antoine-20260915.md).

## Temperature-dependent liquid density

Use `liquid_density_model='polynomial'`, `liquid_density_coefficients=a,b,c,d`
and `liquid_density_temperature_range=Tmin,Tmax` for
`rho(T)=a+T*(b+T*(c+T*d))` in kg/m3, with T in kelvin. This is the pinned
PelePhysics `MPLiqProps::rho_i` convention. The coefficients must already use SI
units. `liquid_density` remains the constant-model input and initialization
reference; the polynomial supplies the actual density when selected.

Injection converts specified diameter to mass at the injection temperature.
Drag, frozen-film evaporation and ETAB/KH-RT use the current droplet density.
Evaporation freezes density over each substep, while the resulting diameter
uses the updated temperature. Heating can therefore change volume without
creating liquid mass. The formation-energy reference remains unchanged; liquid p/rho and surface-energy corrections are
still neglected in this dilute point-particle model.

Coefficients and bounds participate in restart identity. Evaluating outside
the interval or obtaining density outside 1e-3--1e5 kg/m3 rejects the update;
a polynomial is not guaranteed positive over its entire interval by setup.
The default constant-density path retains its checkpoint context.
`cases/spray_3d/water_density.nml` uses an illustrative linear polynomial to
exercise the implementation, not calibrated water properties. See
[focused evidence](validation/liquid-density-20260915.md).

## Temperature-dependent liquid viscosity for breakup

`liquid_viscosity_model='inverse-temperature'` selects
`mu(T)=a+b/T+c/T**2+d/T**3` with `liquid_viscosity_coefficients=a,b,c,d` and
`liquid_viscosity_temperature_range=Tmin,Tmax`. Use kelvin and coefficients
that return Pa s. The convention matches pinned PelePhysics `MPLiqProps::mu_i`.
ETAB damping and KH-RT's Ohnesorge number use this viscosity at the current
droplet temperature. The model requires an enabled breakup model. It can be
combined with temperature-dependent density.

The default remains `liquid_viscosity_model='constant'`, using
`liquid_viscosity`; its new coefficient/range arrays must be zero. With the
inverse-temperature model, the scalar is retained in input identity but the
fit determines viscosity. Evaluations outside the temperature range or outside
1e-10--1e3 Pa s fail transactionally. Fit coefficients and bounds participate
in restart identity. No positivity proof across the whole fit range is implied.

`cases/spray_3d/water_viscosity.nml` exercises combined density/viscosity and
KH-RT using illustrative coefficients, not a calibrated water model. See
[focused evidence](validation/liquid-viscosity-20260915.md).

## Near-boiling vapor-state temperature cap

`boiling_model='surface-cap'` with `critical_temperature=Tc` enables the pinned
PelePhysics vapor-state convention `Ts=min(Tdrop,0.99*min(Tboil(P),Tc))`.
The pressure-dependent boiling temperature comes from the existing constant
reference latent heat and Clausius--Clapeyron equation. The cap affects the
saturation-pressure evaluation; the droplet temperature, energy accounting,
film transport temperature, and density/viscosity fit bounds retain their
existing meanings.

This option currently requires `vapor_pressure_model='clausius'` and Tc>Tref.
It retains the default `boiling_model='reject'` behavior when absent. A
nonpositive inverse boiling temperature is rejected. The critical temperature
and selected model are part of restart identity. The example
`cases/spray_3d/water_surface_cap.nml` starts water droplets at 400 K.

This is a near-boiling numerical vapor-state treatment, not a heat-limited
boiling-rate, supercritical-fluid or bubble-growth model. This application mode
still uses the supplied constant reference latent heat.
See [source and evidence](validation/surface-cap-20260915.md).

The library now separately provides `watson_latent_heat` and
`pelemp_boiling_temperature` in `spray_parcel_mod`. Both return zero and
`ok=.false.` on invalid input. They accept distinct latent-reference and normal
boiling temperatures; the latter is defined at 101325 Pa. Inputs use K, J/kg,
Pa and kg/kmol. The Watson function accepts only subcritical temperatures;
the boiling function caps the inferred temperature at Tc and rejects
nonpositive inverse temperatures. The `pelemp-cap` mode below
connects them to the ordinary EOS vapor-state calculation.
See [kernel checks and limits](validation/watson-properties-20260925.md).

### PeleMP vapor-state cap

Select `boiling_model='pelemp-cap'`, provide `normal_boiling_temperature` at
101325 Pa and `critical_temperature`, and set `reference_temperature` to the
temperature of the supplied `latent_heat`. Both reference temperatures must
be below Tc; `reference_pressure` must be 101325 Pa in this mode. The normal
boiling temperature is independent of the liquid energy reference temperature.
The default `reject` and legacy `surface-cap` modes require the new field to
remain zero.

Watson latent heat evaluated at the normal boiling point supplies Tboil(P).
For Clausius saturation pressure, the latent heat is instead
`hgas(Tdrop) - liquid_energy(Ts)`, where `Ts=min(Tdrop,0.99*Tboil(P))`.
For constant liquid heat capacity this matches the pinned ordinary EOS
PeleMP formula. Polynomial heat capacity uses PeleF's integrated liquid energy.
Antoine evaluation also uses Ts, with the fit bounds applied to Ts.
Stored liquid energy and the conservative transfer update retain their existing
definitions; this is not full reference transfer-integrator parity.

Regular input uses `normal_boiling_temperature`,
`port_normal_boiling_temperatures` or `mixture_normal_boiling_temperatures` with
the corresponding boiling-model and critical-temperature fields. MPI AMR's
`amr_spray_properties` and EB's `eb_spray` groups accept the mixture array.
Unused entries must be zero. All active settings bind checkpoint context.
The expanded AMR/EB namelist serialization may reject checkpoints from older
builds; no cross-build compatibility is implied. Fixed/selected regular
conservation and restart checks are recorded in
[the integration evidence](validation/pelemp-cap-integration-20260925.md).

## Multiple injection ports

Set `injection_port_count=1..16`, then provide `injection_port_centers(3,n)`,
`injection_port_axes(3,n)` and `injection_port_fractions(n)`. In a namelist,
list x/y/z for the first port, then x/y/z for the next. Fractions must be
positive and sum to one (roundoff is normalized). Unused entries must be zero.
Each active center must lie inside the domain and each axis must be nonzero.

`initial_liquid_mass` and `mass_flow_rate` are totals over all ports. Each port
receives its fraction of the event mass and creates `parcel_count` parcels.
IDs advance in port order, preserving deterministic restart and parent counts.
By default the schedule, liquid, temperature, speed, cone angle and diameter
distribution are shared. Independent schedules are available below; distinct
fuels per port and resolved nozzle flow remain outside this path. `maximum_parcels` bounds the combined allocation.

The default count zero uses the original single `injection_center`/`injection_axis`
path and checkpoint context. With active ports, their centers/axes replace the
single geometry; all port data participate in restart identity. See
`cases/spray_3d/water_ports.nml` and [focused evidence](validation/injection-ports-20260915.md).

## Independent port schedules

Set `independent_port_schedules=.true.` with active injection ports, then supply
`port_mass_flow_rates`, `port_injection_starts`, `port_injection_ends` and
`port_injection_intervals` in port order. The global `mass_flow_rate` must be
zero. Each active port needs a positive interval and a nonempty time window.
A zero rate disables scheduled injection at that port. Unused entries are zero.
The aggregate input is limited to one million scheduled events.

Each port injects `rate*min(interval,end-event_time)` and creates `parcel_count`
parcels per event. Events use start+n*interval and occur strictly before end.
The solver processes all due events before stepping and breaks exact timestamp
ties by port number. Initial liquid mass still follows `injection_port_fractions`;
those fractions do not multiply the independent scheduled rates. Material remains shared. Speed, temperature, diameter and cone angle can be
overridden per port as described below.

Restart preserves the aggregate event count and reconstructs the per-port
counts by replaying the deterministic schedule. It checks injected mass,
parent counts, future events and missed events before continuing. Schedule
arrays are included in checkpoint identity. Replay is bounded by the input's
one-million-event limit; large restart performance is not qualified.
See `cases/spray_3d/water_independent_ports.nml` and
[focused evidence](validation/independent-port-schedules-20260915.md).

## Per-port injection conditions

`injection_port_speeds`, `injection_port_diameters`,
`injection_port_temperatures` and `injection_port_cone_angles` override the
corresponding global values for active ports. An entry of -1 (the default)
inherits its global value; unused entries must remain -1. Units are m/s, m, K
and radians. The existing range/material checks apply to each resolved port.

Initial injection and both shared and independent schedules use the resolved
conditions. Different temperatures and diameters produce the corresponding
single-drop mass and multiplicity, while each event retains its specified
mass allocation. Momentum and energy ledgers use each generated parcel's own
velocity and temperature. By default the liquid, distribution shape and truncation bounds
remain shared. For a distributed diameter input, the per-port diameter replaces
the same scale parameter as `drop_diameter`.

Resolved override values participate in restart identity. Without overrides,
the preceding port checkpoint context remains unchanged. See
`cases/spray_3d/water_port_states.nml` and
[focused evidence](validation/injection-port-states-20260915.md).


## Different pure liquids by injection port

`port_vapor_names(16)`, `port_liquid_densities(16)`, `port_liquid_cps(16)`,
`port_reference_temperatures(16)`, `port_reference_pressures(16)` and
`port_latent_heats(16)` override each active port's liquid. Units are SI;
blank names and numeric -1 inherit the global values. Unused slots must keep
those defaults. Every vapor name must exist in the configured gas mechanism;
formation energy is initialized separately from that species' thermodynamics.
Specifying any override enables material IDs equal to the port number.

Initial injection and both schedule modes retain this identity. Droplet mass,
heat transfer, evaporation, gas-species feedback, boundary budgets and system
energy use the selected liquid. Computational splitting copies the identity.
Parcel CSV adds `material_id`; opt-in checkpoint records retain it and reject
changed liquid inputs. Default single-liquid checkpoint/output formats remain.

Each port can also select its density and saturation-pressure models and
near-boiling cap, as described below. Breakup viscosity and surface tension can also be overridden per port;
model selection and the remaining breakup controls stay shared. The low-level material
pool API accepts individual `spray_liquid` records; totals require validated
parcel IDs and materials, as supplied by injection/coupling or checked after
application restart. This is coexisting pure-liquid droplets, not a liquid
mixture with composition-dependent evaporation.

`cases/spray_3d/port_materials.nml` uses H2O and H2O2 in inert nitrogen with
chemistry disabled. Its second-liquid properties are illustrative integration
inputs, not a calibrated peroxide model.
See [focused evidence](validation/port-materials-20260915.md).


### Property models per port

The following arrays select individual liquid models (up to 16 ports):

| Input | Shape | Default |
| --- | --- | --- |
| `port_density_models` | 16 | `constant` |
| `port_density_coefficients` | 4,16 | 0 |
| `port_density_temperature_ranges` | 2,16 | 0 |
| `port_vapor_pressure_models` | 16 | `clausius` |
| `port_antoine_coefficients` | 4,16 | 0 |
| `port_antoine_temperature_ranges` | 2,16 | 0 |
| `port_boiling_models` | 16 | `reject` |
| `port_critical_temperatures` | 16 | 0 |

Select `polynomial`, `antoine` or `surface-cap` as for the corresponding global
model. Coefficients, SI units, bounds and rejection rules are identical.
For example, `port_density_coefficients(:,2)=1150,-.5,0,0` and
`port_density_temperature_ranges(:,2)=280,500` supply the second port's
illustrative linear density fit after selecting its `polynomial` model.

Any nondefault port model input activates per-port material identity even if
all vapor names inherit the global name. Port model fields have their own
explicit defaults; they do not inherit global fitted coefficients or model
selectors. This prevents a fit for one global liquid from being applied to a
different port liquid implicitly. Unused slots retain the defaults. A cap still
requires Clausius pressure and a critical temperature above the liquid reference
temperature; Antoine and cap cannot be combined for the same material.
All port model fields participate in restart identity when enabled. Default
material contexts remain compatible. See
[focused evidence](validation/port-property-models-20260915.md).


### Breakup properties per port

`port_surface_tensions(16)` and `port_liquid_viscosities(16)` override surface
tension (N/m) and dynamic viscosity (Pa s). Numeric -1 inherits the global
value. `port_viscosity_models(16)` defaults to blank, inheriting the complete
global viscosity model and fit. An explicit `constant` or `inverse-temperature`
selector instead uses that port's `port_viscosity_coefficients(4,16)` and
`port_viscosity_temperature_ranges(2,16)`, both defaulting to zero. Coefficients
without an explicit selector are rejected. Inverse-temperature units, validity
bounds and the evaluation formula match the global viscosity model.

These inputs require active ports and ETAB or KH-RT breakup. Keep valid global
breakup settings as inheritance defaults. Unused port fields must retain their
defaults; injection temperature must lie in the chosen viscosity fit's domain.
Later out-of-domain candidates retain the existing whole-step rollback behavior.
Any override enables material identity, including when vapor species are shared.
Both coupling half steps select breakup properties by that identity. Computational
splitting preserves identity; model selection, calibration constants and global
storage/substep controls remain shared in application inputs. Resolved properties
are included in restart identity. With no overrides, preceding formats remain.

The core coupling API accepts an optional `breakup_options` array alongside the
material pool, with one entry per material and the same breakup model as global
options. Only its breakup fields are used; global options still control transfer,
substeps and computational splitting. See
[focused evidence](validation/port-breakup-properties-20260916.md).


## Temperature-dependent liquid heat capacity

Use `liquid_cp_model='polynomial'`, `liquid_cp_coefficients=a,b,c,d` and
`liquid_cp_temperature_range=Tmin,Tmax` for
`cp(T)=a+b*T+c*T**2+d*T**3` in J/(kg K), with kelvin temperatures. The default
`constant` uses `liquid_cp` and requires zero coefficient/range arrays.
The polynomial must remain within [1,1e8] J/(kg K) over its complete interval;
validation checks endpoints and all interior extrema. Tref must lie inside
that interval. Coefficients are bounded by 1e10 in absolute value; temperatures
are within [1,1e5] K. Injection and subsequent live droplets cannot leave the
fit interval. Out-of-range update candidates retain whole-step rollback.

Liquid energy uses the exact polynomial integral from Tref, anchored to the
existing gas enthalpy minus reference latent heat. Both the heating/evaporation
update and its substep estimate use cp at the substep's initial droplet
temperature. This is a frozen-property time approximation; it does not exactly
integrate the nonlinear temperature ODE. Discrete gas feedback uses the change
in integrated liquid energy, preserving all-phase total energy. Clausius
saturation pressure still uses the prescribed reference latent heat; this
change does not add a variable-latent-heat saturation model.

For separate port materials use `port_cp_models(16)`,
`port_cp_coefficients(4,16)` and `port_cp_temperature_ranges(2,16)`. These
independently default to constant/zero, like other per-port property models;
constant cp inherits the existing per-port/global scalar cp input. Any nondefault
port cp field enables material identity. Unused slots retain defaults. All fits
and intervals enter checkpoint identity; default constant-cp formats remain.

The tests use an illustrative `cp=4000+4*T` fit over 300--600 K, with no claim
of measured liquid properties. See [focused evidence](validation/liquid-heat-capacity-20260916.md).


## Ahamed wall-impact correlations

`spray_wall_impact_mod` ports three routines from the pinned PelePhysics
`AhamedSplash.H`: impact-regime selection, five splash-velocity components and
four azimuthal mass fractions. Regimes cover rebound, deposition, dry splash,
thermal breakup, spreading and wet splash. The API accepts the dimensionless
impact parameter `Kv`, wall/boiling temperature ratio, impact angle (30--90
degrees in radians), critical impact parameter and film-thickness/drop-diameter
ratio. Velocity inputs/outputs use m/s. Azimuthal fractions are fractions of
incident mass, summing to the supplied splashed fraction when accepted.

These correlation routines are used by the optional wall-impact coupling below.
Child and film-patch constructors are also exposed separately. Lateral wall-film
transport remains open. Inputs are checked for finite/bounded
values. Negative mass fractions produced by the reference fits return failure
and zero mass outputs; they are not clipped into a different distribution.
In particular, thick-film azimuth fits can produce negative fractions even
within the reference angle range. Velocity outputs preserve signed correlation
values; classification and child generation must interpret them together.

See [pinned source](../references/spray_wall_impact.json) and
[focused C++ comparison](validation/wall-impact-kernel-20260916.md).


### Splash children and deposited film patches

`construct_splash_children` constructs four descendants for a supplied splashed
fraction. Their droplet masses follow the four azimuthal fractions; multiplicity,
temperature and material identity are inherited. IDs start from a supplied unused
ID and `parent_id` retains the incident parcel ID. Deformation and RT timers reset.
Velocities use the reference tangent basis and splash components; children start
one child diameter outward from the impact point. A coordinate-axis fallback
handles the reference's degenerate tangent case. Negative outgoing normal speed
or invalid mass distributions reject without replacing an existing child array.

`construct_wall_film` constructs a stationary cylindrical patch from a supplied
deposited fraction, incident Weber/Reynolds numbers and contact angle in radians.
It retains liquid mass, temperature, material ID, location, normal, film diameter
and thickness. Geometry uses the pinned modified Pasandideh-Fard expression for
the coalesced parcel population. These are individual deposited patches. The
optional thermal/evaporation update below evolves each patch; merging and lateral
film motion remain open.

Both routines return a signed five-component wall budget: mass, three momentum
components and total energy, positive toward the wall. Mass remains entirely in
children or the film (wall-budget mass is zero). Momentum and energy are the
incident selected fraction's inventory minus the constructed inventory. The
caller must remove that selected fraction from the parent exactly once. The
routines leave the input parent unchanged and commit destination state only on
success. The constructors alone do not handle gas feedback, boundary crossings,
domain ownership or coupled restart; the coupled-step and application layers
below supply those connections for regular stationary walls.
See [state/conservation checks](validation/wall-impact-state-20260916.md).


### Complete standalone impact event

`resolve_wall_impact` accepts an incident parent, impact point, outward unit
normal, dimensionless impact parameters, contact angle and child ID range.
It classifies the event, generates children/deposition, updates the parent and
returns the signed wall budget in a single transaction. Rebound retains the
parent with the pinned velocity restitution coefficient; thermal breakup retains
the unsplashed parent with the reference mean splash velocity. Deposition and
spreading put all mass in a film patch. Dry splash divides mass between children
and deposition. Retained parents are placed one remaining-droplet diameter
outward from the supplied impact point. The caller remains responsible for
remaining flight time, domain checks and globally unused child IDs.

For wet splash, the pinned classifier leaves the deposition fraction at zero.
This event API explicitly deposits `1-ms`, where `ms` is its splashed fraction,
so incident mass closes with children and film. It does not entrain mass from
an existing film. This is an explicit closure extending the incomplete reference
branch, not a claim of exact wet-impact parity or physical validation.

Invalid distributions, inward child velocities and invalid film geometry reject
the entire event, preserving parent, child array and film output and zeroing the
wall budget. Successful calls replace the film/children output with this event's
results (empty when absent); a caller must append or accumulate them into its
persistent state. Regular gas coupling and checkpoint persistence are wired into
the application below. Lateral film transport remains a future extension; the
pinned particle mover excludes films from position updates
([source audit](validation/film-reference-scope-20260918.md)). See
[event checks](validation/wall-impact-event-20260916.md).


### Coupled-step wall impact API

`advance_spray_coupling_3d` and `advance_spray_reactive_3d` accept optional
`wall_options` (`spray_wall_options`) and an allocated `wall_films` array.
Enable individual faces with `enabled(6)`, provide wall temperatures (K), a
contact angle (radians), and a maximum film-patch count. Enabled faces must be
slip, no-slip or isothermal gas walls. Coupled impact currently requires stationary
walls. Surface tension and viscosity come from the selected per-material breakup
properties, or global spray options even when aerodynamic breakup is disabled.

The first crossed plane is identified from the parcel substep's position segment.
Impact is applied at that substep's end; retained parents and children are placed
inside the domain and advanced for the remaining coupling interval. This shares
the current first-order particle substep approximation; exact hit-time integration
and simultaneous corner-impact resolution are not claimed. New descendants count
against existing parcel/substep limits. Invalid impact, storage overflow or later
chemistry/flow failure rolls back gas, all parcels, film patches and budgets.

Boiling temperature is inverted from the selected saturation model at sampled
gas pressure: Clausius with reference latent heat, or Antoine strictly within its
fit range. Kv uses the incident normal speed, current liquid density, viscosity
and diameter; Kcrit uses gas pressure converted from Pa to bar. Existing stationary
patches covering the impact point contribute their summed heights to the wet-film
ratio. Patches are appended individually. Their optional thermal/evaporation
update is described below; they are not merged.

`spray_system_totals(...,wall_films=...)` includes film mass, material-resolved
species inventory and integrated liquid energy. The boundary ledger includes
signed wall impulse/energy exchange; deposited liquid remains in system mass.
This capability is also exposed by the application inputs below. See
[focused coupling checks](validation/wall-impact-coupling-20260916.md).


### Wall-impact application and restart

Set `wall_impact_faces(6)` in x-low/x-high/y-low/y-high/z-low/z-high order to
enable Ahamed impacts on chosen stationary wall faces. Specify positive
`surface_tension` and `liquid_viscosity` (or per-port overrides), even with
`breakup_model='none'`. The existing `wall_temperature(6)` supplies impact wall
temperature, including on adiabatic gas walls. `wall_contact_angle` defaults to
pi/3 radians; `maximum_wall_films` defaults to one million patches. Aerodynamic
breakup can remain disabled while wall impacts generate descendants.

Impact-enabled runs write `.films.csv` with material ID, location, normal, mass,
temperature, diameter and height. History adds `film_mass`; `liquid_mass` remains
airborne parcel mass and `system_mass` includes gas, parcels and films. Parcel
CSV adds `parent_id` even without aerodynamic breakup or computational splitting.

The opt-in `_FILM_V1` checkpoint suffix stores the film count, all film records
and parcel lineage. Restart validates material identity, physical bounds,
attachment to enabled walls, mass/volume consistency and aggregate inventories.
Changes to impact settings reject restart. Film-free formats retain their prior
layout. The writer rejects nonempty supplied films unless film serialization is
explicitly enabled, preventing accidental loss of deposition state.

`cases/spray_3d/wall_splash.nml` runs a 300 m/s water parcel into the low-z wall,
producing four descendants and a static deposited patch. Liquid properties are
illustrative and gas chemistry/transport and evaporation are disabled to isolate
impact bookkeeping. See [CLI/restart evidence](validation/wall-impact-application-20260916.md).


### Stationary film heating and evaporation

Set `wall_film_evolution=.true.` and supply `wall_liquid_conductivity` in
W/(m K) for the default constant model. For temperature dependence, set
`wall_conductivity_model='polynomial'`, `wall_conductivity_coefficients(4)` and
`wall_conductivity_temperature_range(2)`. The pinned PelePhysics cubic is
`lambda(T)=a+b*T+c*T**2+d*T**3`, with T in kelvin and coefficients in SI units.
The fit must stay positive and at most 1e4 W/(m K) throughout its stated interval;
endpoints and interior extrema are checked. Film updates outside the interval
are rejected transactionally, including candidate temperatures.

For explicit injection ports, `port_wall_conductivities(16)` overrides the scalar
(default -1 inherits it). `port_wall_conductivity_models(16)` defaults to blank,
which inherits the complete global model. An explicit `constant` or `polynomial`
uses that port's `port_wall_conductivity_coefficients(4,16)` and
`port_wall_conductivity_temperature_ranges(2,16)` (both default zero).
A scalar override of an inherited polynomial requires explicit `constant`.
The global default must remain valid even when ports override it. Unused port
settings are rejected; overrides assign distinct material IDs even if vapor
species and other liquid properties match. The conductivity is evaluated at
film temperature each substep, without mixture averaging or a boiling cap.

For example, a fit can be specified with:
```fortran
 wall_conductivity_model='polynomial',
 wall_conductivity_coefficients=.1,.001,0,0,
 wall_conductivity_temperature_range=300,600,
```
These are illustrative coefficients, not a calibrated liquid dataset. See
[conductivity and restart checks](validation/wall-film-conductivity-20260916.md).

`minimum_film_height` defaults to 1e-7 m (the reference 1e-5 cm cutoff). These
inputs enter checkpoint identity. Existing `heat_enabled` controls gas/wall
heating, and `evaporation_enabled` controls film mass loss as well as droplets.
The default remains static deposition when film evolution is disabled.

For film area A and thickness h, the update uses gas conductance lambda_g/h,
wall heat `2*A*lambda_l*(Twall-Tfilm)/h`, and the pinned Chilton-Colburn mass
transfer correlation with the existing binary-carrier diffusion coefficient.
The gas skin uses the reference coefficient 1/2 for both temperature and vapor
mass-fraction interpolation, with a nonnegative hot-gas temperature difference.
The pinned source comment says one-third, but its actual `rule` constant is 1/2.
Evolving-film checkpoints from the earlier one-third implementation are rejected
by the `wall-film-evolution-v2-half-skin` context; static-film and droplet-only
contexts are unchanged. Local properties/rates are frozen over each
bounded explicit substep. Film heat capacity and liquid energy use the selected
constant or polynomial models. The footprint stays fixed; thickness is recomputed
from remaining mass and current liquid density. An evaporating film below the
cutoff transfers its complete remainder to gas and is removed from the patch list.

Wall heat is a separate signed boundary-energy contribution. Gas receives the
discrete film energy loss plus heat supplied by the wall, and vapor enters the
selected species with zero injected momentum for stationary film. This preserves
the established formation-energy convention; it is not a byte-for-byte port of
the reference energy-source bookkeeping. Gas admissibility and film property
bounds are checked before a substep commits. Newly deposited films only evolve
for their remaining coupling interval; late failure rolls back the whole step.

`cases/spray_3d/wall_film_evaporation.nml` demonstrates deposition followed by
film heating/evaporation. This remains a stationary, fixed-footprint approximation,
without lateral flow, moving-wall coupling, entrainment or calibrated film
qualification. See [focused evidence](validation/wall-film-evolution-20260916.md).


### PelePhysics drag correlation

`drag_model='pelephysics'` selects the pinned `Drag.H` correlation:
Cd is 24/Re for 0 < Re <= 1, 24/Re*(1+Re**(2/3)/6) for
1 < Re <= 1000, and 0.424 above 1000. At zero relative velocity its force
and rate are zero. With `breakup_model='etab'`, Cd is multiplied by
`1+2.6432*deformation`, using the parcel's current deformation. A negative
multiplier rejects the update transactionally. Other breakup models do not
apply this multiplier.

The default `drag_model='schiller-naumann'` retains the preceding correlation.
The selected model enters checkpoint identity, so switching drag models on
restart is rejected. The gas-skin Reynolds number, frozen exponential velocity
update and conservative mass/momentum/energy exchange are shared by both
choices. This ports the reference force law, not its complete time integration
or all spray submodels. See [focused evidence](validation/spray-reference-drag-20260916.md).


### Particle-number-weighted size distribution

Set `diameter_weighting='number'` with positive `diameter_shape` to interpret
the truncated Weibull CDF as a droplet-number distribution, matching the
sampling basis of the pinned PelePhysics Weibull model. Each parcel in a
single injection batch has the same multiplicity (represented droplet count)
and a different total mass. The common multiplicity is the requested batch
mass divided by the sum of its individual droplet masses. This retains the
prescribed total injection mass for any finite parcel count.

The default `diameter_weighting='mass'` retains equal parcel masses and
varying multiplicity. Both choices use deterministic midpoint quantiles and
the existing diameter bounds, rather than random reference sampling. With
shape zero all diameters are identical. `drop_diameter` is the Weibull scale;
it is not generally the arithmetic mean, despite the upstream `mean_dia` name.
The selector applies to initial, scheduled, multi-port and independently
scheduled port batches, and enters checkpoint identity. Switching weighting
on restart is rejected. This does not implement breakup-child generation.
See [focused evidence](validation/spray-number-cdf-20260916.md).


### Normal and lognormal diameters

`diameter_distribution='normal'` or `'lognormal'` enables the corresponding
pinned PelePhysics diameter law. `drop_diameter` is the untruncated arithmetic
mean and `diameter_std` is the untruncated arithmetic standard deviation, both
in meters. For lognormal samples the internal normal parameters are
`width=sqrt(log(1+(std/mean)**2))` and `location=log(mean)-width**2/2`.
Each port's diameter override supplies its own mean; the global standard
deviation and distribution selector provide defaults for the per-port settings below.

The existing `diameter_min` and `diameter_max` truncate and renormalize the
CDF; midpoint quantiles are inverted by bounded bisection. Positive-tail
normal probabilities use the complementary CDF to avoid subtraction from one.
Intervals with numerically unresolved probability are rejected. Positive
standard deviations must lie in [1e-15, 1] m; zero produces the mean diameter
and requires that it lie within the configured bounds. Truncation generally
changes the actual mean and standard deviation of the injected population.

Choose `diameter_weighting='number'` for upstream number-distribution semantics,
or `'mass'` for a mass-CDF interpretation. Both normalize to the requested
injected mass. `diameter_shape` must remain zero for these two laws.
The default `'weibull'` retains the previous shape-zero monodisperse or
positive-shape Weibull behavior and requires `diameter_std=0`.
Selector and standard deviation enter restart identity. See
[normal/lognormal checks](validation/spray-normal-sizes-20260916.md).


### ChiSquared table

`diameter_distribution='chi-squared'` selects the pinned 100-bin `ChiSquared`
distribution. `drop_diameter` (or each port's diameter override) maps to the
upstream nominal `d32` parameter. At bin j=1..100, its cumulative probability
is proportional to `1-exp(-xi)*(1+xi+xi**2/2+xi**3/6)` with `xi=.12*j`;
its diameter is the left edge `.12*(j-1)*d32/3`. This reproduces the source's
discrete table rather than substituting a continuous inverse distribution.
The nominal input does not guarantee the measured Sauter mean diameter after
binning, bounds, weighting and finite sampling.

Bins outside `diameter_min`/`diameter_max` are removed and the remaining
probabilities are normalized. The source's zero-diameter first bin is therefore
excluded by the positive minimum diameter. Empty support is rejected. Midpoint
quantiles select bins; repeated diameters are possible. Select number weighting
for the source sampling basis, or mass weighting for the alternative mass CDF.
`diameter_shape` and `diameter_std` must both be zero. Existing multi-port,
scheduling and restart paths carry the selector. See
[table/restart checks](validation/spray-chi-squared-20260916.md).


### Per-port diameter distributions

Each explicit injection port can override:

- `port_diameter_distributions(16)` and `port_diameter_weightings(16)`;
  blank strings inherit the global selectors.
- `port_diameter_stds(16)`, `port_diameter_shapes(16)`,
  `port_diameter_minima(16)` and `port_diameter_maxima(16)`;
  -1 inherits the corresponding global scalar.

Inheritance is independent for each field. When changing distribution, explicitly
clear an incompatible inherited parameter: for example a normal port under a
positive-shape global Weibull configuration needs `port_diameter_shapes(port)=0`.
Resolved options are validated by the same injector before any run outputs are
created. Nondefault unused-port settings are rejected. Existing port diameters
continue to set the Weibull scale, normal/lognormal mean or nominal ChiSquared d32.

Initial, shared-schedule and independent-schedule batches use these resolved
settings. Each batch conserves its assigned mass fraction or scheduled dose.
Only explicit overrides add the per-port distribution checkpoint context;
changed settings reject restart. The Fortran multi-port API accepts an optional
`spray_size_options` array, one complete record per port. This changes size
statistics, without requiring distinct liquid material IDs.
See [per-port checks](validation/spray-port-sizes-20260916.md).


### Finite circular nozzle outlets

`nozzle_radius` defaults to zero (point injection). A positive value in meters
places injection positions in the circular outlet plane perpendicular to the
specified injection axis. `port_nozzle_radii(16)` defaults to -1, inheriting the
global radius; zero explicitly selects a point outlet for that port.

Following pinned `SprayInjection.H`, the radius fraction is sampled linearly:
`r/R=(i-.5)/N`. This is uniform in radius, not uniform in disk area. Azimuth uses
the existing deterministic ID-based sequence and the same transverse basis as
the velocity cone. Each injected batch retains its previous mass and velocity
rules. The nozzle is a parcel source; it does not carve a gas inflow opening
or replace gas boundary conditions.

The application checks the entire disk's Cartesian extents before running:
`extent_j=R*sqrt(1-axis_unit_j**2)`. The lower edge may touch zero; the upper
edge must remain strictly below the domain extent. This applies even with
periodic gas boundaries. Negative/nonfinite radii and unused port overrides
are rejected. Initial, shared and independently scheduled injection all use
the resolved radii. Positive/overridden radii enter restart identity; changed
geometry rejects restart, while default point contexts remain unchanged.

The Fortran injectors accept optional `nozzle_radius` and, for multi-port use,
`port_nozzle_radii` (complete values, no inheritance sentinel). This ports
finite outlet geometry and radial sampling. Hollow and swirl options are below;
time-varying geometry and reference random trajectories remain open. See
[geometry/restart checks](validation/spray-nozzle-20260916.md).


### Solid/hollow nozzle profiles and swirl

`nozzle_profile='legacy'` retains the original equal-solid-angle velocity cone.
`'solid'` selects the source's linear relation between radius fraction q and
velocity inclination: r=R*q and theta=`cone_angle`*q. `'hollow'` places all
parcels at r=R and uses theta=`cone_angle`+`hollow_spread`*(q-.5). Both new
profiles require positive nozzle radius. Here q=(i-.5)/N, a deterministic
midpoint replacement for the source's random samples. Hollow inclination
spread defaults to zero and is allowed in [0,pi] radians; it can produce
negative or backward inclinations for sufficiently wide spreads, as in the
source. Non-hollow profiles require zero `hollow_spread`.

`swirl_angle` rotates velocity azimuth relative to the outlet position azimuth
without changing speed or axial velocity. It defaults to zero and must lie in
[-pi,pi] radians. All these application angles are radians; the upstream jet
input uses degrees. `cone_angle` remains a half-angle.

Per-port overrides are `port_nozzle_profiles(16)` (blank inherits),
`port_swirl_angles(16)` (-999 inherits), and `port_hollow_spreads(16)`
(-1 inherits). Inheritance is independent: a solid override under a hollow
global profile must explicitly clear any inherited spread. The selectors and
angles enter restart identity. Existing default contexts are unchanged.
The Fortran injectors accept optional `spray_nozzle_model` records, globally
or as a complete array per port. See
[profile/swirl checks](validation/spray-nozzle-profile-20260916.md).


### Time-dependent rate of injection

Set `injection_profile_points` to 2..64, with matching
`injection_profile_times(64)` in seconds and `injection_profile_rates(64)` in
kg/s. Times must strictly increase, cover the entire injection window and lie
in [0,1000] s; rates must lie in [0,1] kg/s. Entries after the declared count
must be zero. `mass_flow_rate` must be zero when the table is supplied.

The injector linearly interpolates between knots and integrates each scheduled
window exactly by trapezoids, including partial knot intervals. A batch is
still deposited at the beginning of its injection window, as in the existing
schedule; this is not continuous-in-time particle creation. Zero-dose windows
advance the schedule without generating parcels. Shared ports split each dose
by their configured fractions and retain their own liquid/size/nozzle settings.
The global table uses the shared schedule. Independent-port tables are configured
separately below; a global table cannot be combined with independent mode.

Table data enter checkpoint identity. Restart reconstructs both accumulated
mass and the number of nonzero batches, including pauses, before accepting
particle/event counts. Changing a knot or rate rejects restart. A zero point
count preserves the previous constant-rate behavior. See
[rate-profile checks](validation/spray-injection-profile-20260916.md).


### Independent port rate histories

With `independent_port_schedules=.true.`, each port may use
`port_profile_points(16)` and corresponding columns of
`port_profile_times(64,16)` / `port_profile_rates(64,16)`.
A zero point count retains that port's constant `port_mass_flow_rates` value.
A positive count (2..64) requires the corresponding constant rate to be zero.
Table and constant-rate ports may coexist.

Each table must cover its own `port_injection_starts`/`port_injection_ends`
window, with increasing times and nonnegative rates in the same ranges as the
global table. Unused points and ports must have zero table entries. Each event
uses the exact integral over its own port interval; zero-dose events advance
the schedule without creating parcels. Simultaneous events retain port-index
ordering. The global table and these independent tables are mutually exclusive.

Tables enter restart identity. Restart replays the merged port event order,
reconstructs individual window doses and checks aggregate injected mass and
nonzero-batch parcel counts. See
[independent profile checks](validation/spray-port-rate-profiles-20260916.md).


### External rate files

The regular serial application's shared or independent tables can alternatively
be loaded from an optional group after `&spray`:

```fortran
&spray_rate_files
 profile_file='shared.dat'
/
```

For independent schedules, use `profile_files(p)='port.dat'` instead of the
shared `profile_file`. The shared and independent modes remain mutually
exclusive. Leave the corresponding inline count/arrays at zero. Each file has
2--64 whitespace-separated time/rate pairs; blank lines and `#`/`!` comments
are allowed. Existing 0--1000 s, 0--1 kg/s and event-coverage bounds still apply.
Unused port paths and inline/file conflicts reject. Relative paths use the
working directory. Files load at startup, including restart; loaded numerical
tables, not paths/comments, bind checkpoint identity. Existing event scheduling
is unchanged. See [serial file checks](validation/spray-serial-rate-files-20260925.md).

### Velocity from injection rate

`velocity_from_rate=.true.` derives scheduled injection speed from
`U=mdot/(rho*pi*nozzle_radius**2*discharge_coefficient)`, following the pinned
`SprayJet::readROI` relation. `discharge_coefficient` defaults to 1 and must be
in [1e-6,1]; `port_discharge_coefficients(16)` defaults to -1, inheriting the
global value. Each active configuration requires a nozzle radius of at least
1e-12 m. Density is evaluated from that port's liquid model at its injection
temperature. Derived speed is checked against the injector's existing limit.

Here mdot is the integrated dose divided by the batch window duration. Shared
ports use their assigned fractions; independent ports use their own constant
or tabulated rates and intervals. The flow-derived speed replaces prescribed
`injection_speed`/port speeds in this mode; cone, hollow and swirl settings still
set its direction. The application probes peak-rate speeds before running.

This mode requires `initial_liquid_mass=0` because an initial inventory has no
associated injection window. Coefficients are otherwise unused and cannot be
set away from their defaults when the mode is disabled. Mode and coefficients
enter restart identity. The default prescribed-speed behavior retains its
context. Each batch still has one speed and is emitted at the window start;
continuous momentum-flux equivalence is not claimed. See
[flow-derived velocity checks](validation/spray-rate-velocity-20260916.md).

## Cartesian particle partitioning

`spray_partition_mod` supplies deterministic ownership and send-buffer ordering for
the forthcoming particle MPI path. Each Cartesian axis assigns remainder cells to
the lowest coordinates. Rank numbering is x-fastest. Physical coordinates must be
inside the half-open box; periodic coordinates wrap for ownership without changing
the stored particle position. Partitions cannot exceed cell counts or one million
ranks in total.

`partition_spray_parcels` accepts globally increasing particle IDs and groups records
by destination rank, preserving ID order within each group. Counts and zero-based
offsets identify contiguous send segments. The complete parcel record and optional
per-drop component matrix move together. Invalid input leaves every caller array
unchanged. MPI communication and receive-side ID ordering are available below; distributed
gas feedback remains open. See [partition evidence](validation/spray-partition-20260916.md).

### MPI particle migration

With `PELEF_ENABLE_MPI=ON`, `migrate_spray_parcels` uses `mpi_f08` all-to-all
communication to move records to Cartesian owner ranks. IDs and lineage use integer
64-bit buffers; the 14 real parcel fields and optional component masses use separate
double-precision buffers. Received records are sorted by ID before committing. All
ranks must supply identical geometry, decomposition and component-column count.
The caller supplies globally unique IDs; duplicates arriving at one rank reject.

Invalid local input, mismatched configuration, receive-capacity overflow or invalid
received records reject collectively before any caller arrays are replaced. This
guarantee concerns application validation, not failed processes or MPI transport
errors. Empty sources and destinations are supported. This migration API does not
yet distribute gas coupling, film patches, injection IDs or application checkpoints.
See [MPI evidence](validation/spray-mpi-migration-20260916.md).

### Collective particle ID reservation

`reserve_spray_ids` reserves a contiguous interval for each rank from a shared
64-bit `last_id`. Rank order determines interval order. Each rank may request zero
to one million IDs; a zero request returns `first_id=0`. On success every rank
receives the same updated high-water mark. Negative requests, mismatched existing
marks and overflow reject collectively without changing either output value.

The caller must preserve the high-water mark in application restart state and
reserve before constructing injected or split particles. The API does not infer
that mark from locally owned records. Assignment is deterministic for a given rank
layout and request vector; changed-rank particle identity parity is not claimed.
Injection/breakup driver wiring remains open. See
[ID evidence](validation/spray-mpi-ids-20260916.md).

### Collective computational splitting

`split_mpi_spray_parcels` applies the existing multiplicity cap on local candidates,
reserves globally disjoint child IDs and commits particle/component state and the
shared ID mark together. Children inherit per-drop component masses and parent IDs;
KH initial count and shed volume are partitioned with multiplicity. If any rank
cannot split within capacity, no rank commits or consumes IDs. Input IDs must be
locally increasing and no larger than the shared mark. The component table is a
caller-owned common convention; column counts must agree across ranks.

This is computational splitting of represented drop counts, not aerodynamic
breakup. MPI injection and distributed gas/film coupling remain to be connected.
See [split evidence](validation/spray-mpi-split-20260916.md).

### Collective cone injection

`inject_mpi_cone_parcels` combines global ID reservation, the existing cone injector
and a collective commit. Each rank supplies its own injected mass, count, position,
velocity and size/nozzle options. Zero-count ranks supply zero mass. Optional common
component tables and rank-local fractions create mixed droplets using density and
energy evaluated at injection temperature. Component columns append with particles.
Returned `injected_totals` are local mass, momentum, energy and gas-species amounts;
the caller sums or records them as appropriate.

A failed local geometry/property check or capacity check preserves every rank's
particles, component matrix and shared ID mark, returning zero injection budgets.
Successful injection appends increasing global IDs but does not automatically move
particles to spatial owners; call migration after injection. Distributed scheduling
and the MPI CFD application remain open. See
[injection evidence](validation/spray-mpi-injection-20260916.md).

### Distributed gas-field sampling

`sample_mpi_spray_field` takes each rank's Cartesian cell block and local particle
positions. It sends the eight interpolation-cell indices to their owning ranks and
returns only the requested field values, then evaluates the ordinary particle
stencil locally. No root or replicated whole-grid field is required.

All ranks must agree on global cells, spacing, lower corner, Cartesian parts,
boundary faces and field width. `spray_cartesian_block` supplies zero-based first
cells and local extents. Empty requesting ranks are supported. This bounded API
accepts up to 125000 local positions, 256 fields and one million incoming cell
requests per rank. Invalid configuration, positions, field values or received
indices reject collectively while preserving previous sample arrays. These are
interpolated conserved fields; thermodynamic recovery remains the caller's job.
Distributed feedback deposition and time stepping remain open. See
[sampling evidence](validation/spray-mpi-sampling-20260916.md).

### Distributed feedback deposition

`deposit_mpi_spray_field` distributes local particle inventory increments through
the same eight-cell stencil. It routes indices and weighted contributions to cell
owners, divides by cell volume and adds to local field candidates. It needs no root
or whole-grid copy. Sources are accumulated in sending-rank order, then local
particle/stencil order; changed-rank floating-point parity is not implied.

The geometry, field-width and request limits match sampling. Empty source ranks
are supported. Invalid input, request routing or floating-point overflow rejects
collectively and preserves every local field. This is a generic conservative
field operation: gas positivity, species closure and temperature recovery must be
checked by the future coupled driver before committing its particle candidates.
See [deposition evidence](validation/spray-mpi-deposition-20260916.md).

### Thermodynamically checked distributed feedback

`deposit_mpi_spray_gas` applies inventory increments to distributed conserved gas
blocks and recovers temperatures with the existing reactive thermodynamics. It
validates the initial state, deposits into a private candidate, and accepts both
state and temperature only when every owner recovers successfully. Invalid density,
species or thermal energy on one owner preserves both arrays on all ranks. The
ordered species table must be identical on every rank. Geometry and request limits
are inherited from field deposition.

The caller still owns particle candidates, time-step retries and the combined
particle/gas commit. This API does not itself advance particles or chemistry.
See [gas feedback evidence](validation/spray-mpi-gas-20260916.md).

### Distributed periodic mixed-particle coupling

`advance_mpi_mixture_coupling` connects owner-routed sampling, mixture drag,
heating/evaporation, gas deposition and temperature recovery over a supplied
interval. Ranks choose a common minimum of particle drag, thermal, component
loss and flight limits. Each trial advances all local particles from the same
sampled gas, then accumulates feedback collectively. A rejected trial halves the
common step (up to 30 attempts). Exhaustion or final migration failure preserves
all caller gas, temperatures, particles and component inventories; `substeps` is
zero on failure. Accepted intervals remove evaporated particles and migrate
survivors, retaining IDs and component columns.

This is a synchronous particle update, so equality to the serial coupler's
particle-by-particle feedback order is not claimed. Each rank needs identical
ordered species/transport/component tables and interval; it can start with no
particles. Existing sampling limits apply (125000 local particles, at most 255
conserved fields because temperature occupies one extra sampled field). Only
periodic boundaries and `breakup_model='none'` are currently supported; unsupported
breakup is rejected. Wall events, distributed breakup, gas hydro/chemistry stepping,
injection scheduling and the MPI CLI remain to be connected. See
[coupling evidence](validation/spray-mpi-coupling-20260916.md).

### Distributed particle and reaction splitting

`advance_mpi_mixture_reaction` performs a half-interval periodic mixed-particle
update, cell-local gas chemistry, then the remaining particle half. All stages
operate on private gas, temperature, particles and composition; collective
chemistry failure or a failed second half discards the complete interval. Accepted
particle substep counts are returned only on success. Interval, tolerances and
integrator selection must agree across ranks. Ordered thermodynamic, reaction,
transport and component tables are a caller-owned common convention.

The wrapper uses `advance_spray_chemistry_3d`, accepting native explicit/implicit
chemistry and CVODE when compiled in. Focused MPI evidence currently covers native
implicit chemistry. Without hydro options this is particle/reaction operator splitting. SGS, walls
and breakup are not advanced by this API. The
nonlinear particle substep remains first order despite symmetric ordering. See
[reaction coupling evidence](validation/spray-mpi-reaction-20260916.md).

### Distributed periodic gas advection

`advance_mpi_spray_hydro` exchanges two x-cell layers with adjacent ranks and
uses the existing characteristic PLM fluxes and SSPRK2 update. Decomposition is
`[nranks,1,1]` with at least two x cells per rank and three global cells in y/z.
All directions are periodic. Riemann solver (`rusanov`, `hllc`, `pelec`), limiter
(`minmod`, `mc`), geometry and time step must agree collectively; species tables
are the caller's common convention. The caller chooses a stable time step.
Intermediate stages stay private until all ranks recover valid temperatures.

Supplying both `hydro_solver` and `hydro_limiter` to
`advance_mpi_mixture_reaction` enables particle-half / chemistry-half / hydro /
chemistry-half / particle-half ordering. Gas and particle state commit only after
all stages and final migration succeed. Omitting both retains reaction-only gas
stepping. The distributed application remains open. See
[advection evidence](validation/spray-mpi-hydro-20260917.md).

### Distributed molecular transport

`advance_mpi_spray_transport` advances periodic x-slab gas viscosity, heat
conduction, species diffusion and barodiffusion using the existing ghosted
transport fluxes and SSPRK2. Its four `enabled` entries select these processes in
that order. One gas layer and cell species-outflow limiter factors are exchanged
with adjacent ranks each stage, so each shared face uses the same limited flux.
The returned `minimum_theta` is the minimum accepted factor across ranks/stages;
failure leaves gas/temperature intact and returns its default value of one.
Geometry, time step and flags must agree; species and transport tables remain a
caller-owned common convention. The caller chooses a stable time step.

Optional `transport_flags(4)` in `advance_mpi_mixture_reaction` requires hydro
options when any flag is true and adds transport half-steps around hydro:
P/2-C/2-T/2-H-T/2-C/2-P/2. All stages retain the complete interval transaction.
Gas SGS, walls and breakup remain open. See
[transport evidence](validation/spray-mpi-transport-20260917.md).

### MPI executable

With spray and MPI enabled, `pelef_mpi_spray_3d` reads a root-broadcast `mpi_spray`
namelist, initializes the fixed H2/O2 gas and a configurable two-component cone
population, and advances the coupled periodic x-slab solver. Global hydro/transport
limits choose the time step; rejected complete updates are retried by halving.
History records global budgets and particle count, while rank-local gas/particle
CSVs retain final fields. Existing output files are refused. See the
[case and input guide](../cases/mpi_spray_3d/README.md) and
[application evidence](validation/spray-mpi-application-20260917.md).

The MPI executables currently support timed repeated injection and periodic
flow. Its CSVs do not provide restart; the checkpoint API below does. Independent multiport/rate-profile
injection, SGS and distributed wall/breakup events remain open.

### Checkpoint ID reservation state

`write_spray_checkpoint` and `read_spray_checkpoint` accept optional 64-bit
`last_id`. Supplying it adds `_ID_V1` to the format identity and stores an integer
record after dimensions/component count, before time metadata. The reader must
explicitly request this mode. The watermark must be nonnegative and at least every
live particle ID; it may exceed all live IDs because particles can evaporate.
Overflow, malformed records and smaller watermarks reject without changing state
or the caller's watermark. Existing checkpoint modes remain unchanged when this
argument is omitted. This is the persistence API needed by collective MPI ID
reservation and is used by MPI field/particle collection and CLI restart below. See
[ID checkpoint evidence](validation/spray-id-checkpoint-20260917.md).

### MPI checkpoint collection and redistribution

`mpi_spray_checkpoint` in `spray_mpi_checkpoint_mod` writes or reads the existing
mixed-particle checkpoint with `_ID_V1`. It supports periodic x-slab gas blocks,
mixed particle columns, time/step/event metadata, reference/injected inventories
and the shared ID high-water mark. Root owns file I/O and the temporary whole-grid
buffer; other ranks retain local blocks. Writes gather complete particle records
and sort globally by ID. Reads scatter gas into the current partition and migrate
particles to spatial owners, committing caller state only after all operations
succeed. The rank count is not encoded into the file identity.

`migrate_spray_parcels` now accepts optional `gather_rank`: it sends all complete
records to that rank while retaining geometry and payload checks. Omitting it
retains ordinary spatial migration. The checkpoint routine uses private copies.

The checkpoint API currently handles the no-breakup mixed-particle mode used by
the MPI executable. Its caller must bind geometry, ordered property/mechanism
and evolution identity in `context`; root path/context govern file access.
The temporary root gas buffer is bounded to 100 million real values and particle
collection to one million records. These are bounded root-I/O APIs, not scalable
parallel checkpoint I/O. The MPI CLI now uses these APIs for save/restart. See
[collection evidence](validation/spray-mpi-checkpoint-20260917.md).

### MPI executable restart

`checkpoint_path` saves the last accepted state; `restart_path` restores it before
initial history output and skips initial injection. A model/evolution context binds
geometry, mechanism, ordered liquid properties, injection definition and numerical
controls. Output paths, stop time, maximum-step guard and rank count are excluded
so runs can continue in a new decomposition. Initial reference and injection
budgets remain those of the saved run. See the
[restart guide](../cases/mpi_spray_3d/README.md#save-and-resume) and
[application evidence](validation/spray-mpi-restart-application-20260917.md).

### MPI timed injection

`injection_events`, `injection_start` and `injection_interval` schedule repeated
cone events using the existing `number` and `total_mass` per event. Default values
retain one event at time zero. The step ends at the next pending event and the
application processes events at its final time before saving. Event candidates,
component columns and global IDs commit after injection and migration succeed.
History budgets subtract accumulated injection from gas-plus-particle totals.

Restart preserves the event counter and resumes only future events, including
when the checkpoint is written exactly at an event time. Schedule parameters
enter checkpoint context. Context v2 records initial-gas reference budgets and
separate injected budgets, so earlier v1 application contexts are intentionally
incompatible. Independent ports and time-varying rate profiles remain open for
this MPI executable. See [schedule evidence](validation/spray-mpi-schedule-20260917.md).

### Selected-mechanism MPI application

`pelef_mpi_spray_3d_selected` is generated when MPI, spray and a mechanism bundle
are enabled. Both MPI frontends call `run_mpi_spray_application`; gas arrays,
flux budgets, CSV columns and initial mass fractions follow the loaded species
count. Default H2/O2/N2 fractions are assigned by species name; explicit numeric
`gas_fractions` uses bundle order. Checkpoint context contains the selected bundle
SHA and retains schedule/restart identity checks.

The 21-species FFCM1_Red water/methanol example passes the application, schedule
and same/changed-rank restart checker. Its strict nonnegative-species check exposed
a chemistry closure residual on absent HE; the separate
[roundoff correction](validation/reactive-zero-species-closure-20260917.md) preserves
mass by assigning a tiny negative closure residual to the largest component.
Build/input instructions are in the [MPI case guide](../cases/mpi_spray_3d/README.md),
with exact configuration and evidence [here](validation/spray-selected-mpi-20260917.md).

### MPI mass-flow profiles

The MPI namelist accepts up to 64 ordered time/rate points and integrates each
scheduled event interval with `integrate_spray_injection_profile`. The serial
application now calls the same interval integrator, preserving its trapezoidal
formula. Profiles cover all event intervals; unused entries stay zero. Active
profiles replace fixed event mass and are included in restart identity. Zero-dose
events advance the schedule without consuming IDs. Independent MPI port schedules
remain open. See [profile evidence](validation/spray-mpi-profile-20260917.md).

### MPI EB application impacts and films

The fixed and selected MPI EB spray executables now accept `wall_impacts` and
`walls` in `/eb_spray/`. They connect impact-enabled reactive stepping, stationary
film evolution, composite budgets, per-rank film CSVs, regrid/rebalance migration
and generation-bound V2 checkpoint payloads. Restart retains splash-child IDs
and retired parents through a shared watermark and validates wall settings in
the normalized context. Inadmissible empirical impact output fails immediately.
See [input/output controls](../cases/mpi_eb_spray_3d/README.md#eb-wall-impacts-and-stationary-films)
and [focused application evidence](validation/spray-mpi-eb-impact-application-20260918.md).

MPI EB ports also accept `port_radii`, `port_nozzles` and `port_sizes`, reusing
the existing finite solid/hollow outlet, swirl and bounded size-distribution
kernel. Their normalized values enter checkpoint identity. Point/monodisperse
defaults remain available. See the
[case controls](../cases/mpi_eb_spray_3d/README.md#finite-outlets-and-size-distributions)
and [MPI nozzle evidence](validation/spray-mpi-eb-nozzle-20260918.md).

MPI EB ports can select `port_velocity_from_rate` and per-port discharge
coefficients. Each positive event uses its scaled integrated mass divided by
interval duration, mixture density, outlet area and coefficient to determine
speed. Particle velocities and injected momentum/energy use that same value;
zero-mass events consume no IDs. See
[rate-velocity controls](../cases/mpi_eb_spray_3d/README.md#speed-from-the-port-mass-flow-profile)
and [evidence](validation/spray-mpi-eb-rate-velocity-20260918.md).

MPI EB component tables now expose the existing density/heat-capacity polynomial,
Antoine pressure and Clausius surface-cap controls through `mixture_*` fields.
Port composition and temperature supply mixture properties to injection and
rate-derived speed; source/film operators use the full component records.
Temporary injection inventory liquids retain evaluated constant mixture values,
without reapplying a component's polynomial. See
[component controls](../cases/mpi_eb_spray_3d/README.md#component-property-models)
and [focused evidence](validation/spray-mpi-eb-properties-20260918.md).

MPI specular outflow and empty-population rank-changing restart are covered by
[this application record](validation/spray-mpi-eb-specular-application-20260918.md).

## Serial EB external acceleration

`advance_eb_spray_reactive_allocated` accepts optional `acceleration(3)` and
`source_change(nvar)`. Nonzero acceleration requires the output ledger; each
component must be finite and at most 1e6 m/s² in magnitude. Half kicks bracket
the existing coupled step, changing gas and live-particle momentum and kinetic
energy. Gas contributions use fluid volume fraction times cell volume; covered
cells and stationary wall films receive no kick. Temperature recovery checks
that gas internal energy is unchanged. The returned system impulse/work is
separate from exterior, wall and chemistry inventories, with zero mass/species
entries. All state and ledgers publish only on full-step success. The fixed-array
wrapper also forwards `acceleration` and `source_change`; see
[wrapper evidence](validation/spray-eb-force-fixed-20260918.md). The serial application
accepts `acceleration` in `eb_spray`, prints `Force inventory`, and persists it
separately from wall exchange through restart. See
[application evidence](validation/spray-eb-acceleration-application-20260918.md).
See [focused evidence](validation/spray-eb-acceleration-20260918.md).

## MPI EB external acceleration

`advance_mpi_eb_spray_reactive` now accepts `acceleration(3)` and
`source_change(nvar)`. Acceleration must agree across ranks and be finite with
each component bounded by 1e6 m/s²; nonzero acceleration requires source output.
Half kicks bracket the coupled update. Owner-local fluid cells and live parcels
receive impulse/work; covered solid cells and attached films do not. Coarse
cells under children are updated for consistent donors but omitted from the
composite source ledger. Local inventories are summed once across MPI ranks.
Gas temperature recovery checks unchanged internal energy. Failure on any rank
preserves caller state and ledgers. The MPI application now exposes acceleration
and saves cumulative force separately in the gas checkpoint prefix; its budget
CSV includes a `force` column. See
[application evidence](validation/spray-mpi-eb-acceleration-application-20260918.md).
See [focused evidence](validation/spray-mpi-eb-acceleration-20260918.md).

### Additional LES stress kernels

`les_smagorinsky_3d_mod` now exports `wale_stress` and `vreman_stress` with the
same gradient/density/spacing/coefficient/output layout as `smagorinsky_stress`.
These implement the pinned deviatoric WALE and Vreman algebra; zero gradients
return zero safely. See
[reference and focused algebra checks](validation/les-algebraic-closures-20260918.md).

The regular LES RHS, supplied-halo RHS and SSPRK2 library APIs now accept an
optional trailing `model='smagorinsky'|'wale'|'vreman'`. Omission retains
Smagorinsky. The coefficient argument `cs` means Cw for WALE and Cs otherwise.
The selected viscosity drives stress, heat/species fluxes and the diffusion
stability bound. Invalid names reject even when coefficient or dt is zero.
See
[grid integration checks](validation/les-alternative-grid-20260918.md).

`advance_mpi_spray_les` also accepts the optional `model` argument. It checks
model validity and agreement across ranks before halo exchange and forwards
the model through both RK2 stages. Two/four-rank tests cover periodic and moving
isothermal-wall parity with serial; see
[MPI LES model checks](validation/les-mpi-models-20260918.md).

The serial `advance_spray_reactive_3d` and MPI `advance_mpi_mixture_reaction`
coupled library APIs accept optional `les_model` and forward it through both
LES half steps. The MPI wrapper includes the model in collective control
agreement; omission keeps Smagorinsky. Serial reacting/evaporating conservation
checks cover all three models; see
[coupled LES evidence](validation/les-coupled-models-20260918.md).

### Application LES selection

The regular serial `spray` and MPI `mpi_spray` namelists accept
`les_model='smagorinsky'|'wale'|'vreman'`. Serial `smagorinsky_constant` and MPI
`les_cs` retain their input names and bounds; they mean Cw for WALE and Cs
otherwise. Their zero defaults disable SGS. Model validity is checked even with
zero coefficient. The selected model supplies both SGS half steps and their
existing diffusion rejection bound; rejected coupled steps retry at half dt.
Nondefault models enter checkpoint identity, so changing or omitting one on
restart rejects before history output. Omitted/default Smagorinsky preserves
the previous context. Fixed and selected frontends share these drivers.
See [focused application evidence](validation/les-application-models-20260918.md).

### Dynamic Smagorinsky building blocks

`les_dynamic_smagorinsky_mod` provides `dynamic_smagorinsky_quantities` from
supplied cell-centered velocity/temperature gradients and
`dynamic_smagorinsky_coefficients` from supplied filtered tensors. The latter
returns `Cs^2`, `CI`, `Cs^2/PrT` with pinned regularization and clipping. These
APIs do not perform spatial filtering or add a `dynamic` application choice;
filtering, smoothing and dynamic flux integration remain to be connected.
See [reference, numerical boundaries and algebra checks](validation/les-dynamic-algebra-20260918.md).

`les_filter_3d_mod` now provides `les_filter_weights` for standard box, normalized
Gaussian and identity kernels, and `filter_les_field_3d` for periodic fields or
caller-supplied halos. The box ratio-six kernel is available for dynamic coefficient
smoothing. Output is preserved on rejection. These building blocks still require
connection to dynamic LES field evaluation and fluxes; see
[filter evidence and reference scope](validation/les-filter-20260918.md).

`compute_dynamic_les_coefficients_3d` now connects these building blocks on a
periodic regular grid: conserved-state filtering, EOS recovery, centered
primitive gradients, grid/test contractions and ratio-six box smoothing.
It returns the three dynamic coefficient fields transactionally; dynamic face
fluxes, time integration, application selection and nonperiodic/MPI/AMR field
paths remain to be connected. See
[coefficient-field checks](validation/les-dynamic-field-20260918.md).

`compute_dynamic_les_rhs_3d` now forms periodic dynamic SGS momentum/energy flux
divergence, including isotropic stress and dynamic thermal transport. It uses
arithmetic face interpolation and preserves RHS on rejection; mass/species SGS
fluxes are zero as in the pinned dynamic model. Periodic serial/MPI stepping and
application selection are connected below; nonperiodic and AMR paths remain open. See
[dynamic spatial flux checks](validation/les-dynamic-flux-20260918.md).

The periodic serial RHS also accepts optional `fx`, `fy`, `fz` arrays, supplied
together with shapes `(nv,nx+1,ny,nz)`, `(nv,nx,ny+1,nz)`, `(nv,nx,ny,nz+1)`.
Index 1 is the low face; periodic end faces are duplicated. These are flux
densities, not time/area integrals. Mass/species fluxes are zero. RHS and faces
publish only after successful evaluation; invalid requests preserve them.
The cropped `active_x` MPI workspace mode cannot request full-domain faces.
This exposes existing periodic fluxes for future conservative coupling; it does
not yet supply AMR halos, reflux or Dynamic LES AMR stepping. See
[face-output checks](validation/les-dynamic-face-output-20260925.md).

`compute_dynamic_les_rhs_ghosted_3d` now evaluates one rectangular patch from
caller-supplied conserved state/temperature with `ng` donor layers on all six
sides. Supply donor values at the patch spacing. For filter radius `r`, require
`ng >= r+5` (one derivative layer, three smoothing layers, one face-neighbor
layer); `ng` is bounded to 128. RHS dimensions define the interior; state and
temperature dimensions must be the interior plus `2*ng` on each axis. All three
face outputs use the shapes above and publish with RHS only on success.
The implementation evaluates a padded periodic workspace and crops its interior;
the checked halo keeps its artificial seam outside the requested stencil.
Every workspace cell/filter state must remain thermodynamically admissible.
This first implementation overcomputes the outer workspace. Optional
`maximum_diffusivity` bounds the faces that contribute to the returned interior
divergence, including its physical boundary faces, publishing only on success.
The unused outer workspace is excluded from this timestep estimate
([startup correction](validation/dynamic-active-face-timestep-20260925.md)). The caller still owns physical/coarse-fine
donor construction, flux synchronization/reflux and time stepping. See
[ghosted-patch checks](validation/les-dynamic-ghosted-20260925.md).

`advance_dynamic_les_ssprk2_3d` now advances the periodic dynamic operator with
stage-wise coefficient recomputation, EOS recovery and transactional publication.
Each stage applies a frozen-coefficient diffusion estimate with a 0.2 diffusion
number gate, including isotropic stress and thermal transport. This estimate is
not a general nonlinear stability proof. Regular physical boundaries and their
SSPRK2 budgets are now connected; see
[physical operator checks](validation/regular-dynamic-boundaries-20260925.md) and
[time-integration evidence and initial-test diagnosis](validation/les-dynamic-step-20260918.md).

The serial coupled `advance_spray_reactive_3d` API accepts `les_model='dynamic'`
with periodic or supported physical boundaries. Optional `les_test_filter`, `les_test_ratio` and
`les_grid_ratio` default to box, two and one. Dynamic SGS runs in both half steps;
`cs=0` disables it and positive cs enables it without scaling the computed
coefficients. Fixed prt/sct do not set dynamic heat/species transport. Invalid
filters reject even when cs is zero. Both LES half steps return their boundary
exchange into the coupled ledger, with whole-step rollback. See
[physical application checks](validation/regular-dynamic-spray-20260925.md) and
[coupled dynamic checks](validation/les-dynamic-coupled-20260918.md).

### Serial dynamic LES input and restart

Both serial spray executables now accept these `spray` namelist entries:

```fortran
 les_model='dynamic', smagorinsky_constant=.1,
 les_test_filter='box', les_test_ratio=2, les_grid_ratio=1,
```

Periodic, reservoir-inflow, outflow/pressure-target, slip-wall and
adiabatic/isothermal moving-wall boundaries are supported. Positive `smagorinsky_constant` enables
this model; its magnitude and fixed turbulent Prandtl/Schmidt values do not set
the dynamic coefficients. Zero disables SGS. Box/Gaussian/identity test kernels
use the filter API's valid ratios; the grid ratio is 1–64. Nondefault filter
controls require the dynamic model. Model and filter settings are checked on
restart, including settings omitted back to their defaults. MPI dynamic input
remains unsupported. See [fixed/selected application checks](validation/les-dynamic-application-20260918.md).

### Distributed LES filter

`les_mpi_filter_mod::filter_mpi_les_field_3d` applies supplied normalized weights
on periodic x slabs. One-plane neighbor propagation supports stencils wider than
local slabs without globally replicating fields. Weight/configuration agreement
and rejection are collective; output publishes only after all ranks accept.
Two/four-rank unequal-slab checks match serial filtering exactly. Dynamic MPI
coefficient evaluation and time stepping still require connection; see
[distributed filter evidence](validation/les-mpi-filter-20260918.md).

`les_mpi_dynamic_mod::compute_mpi_dynamic_les_coefficients_3d` now evaluates and
smooths dynamic coefficient fields on periodic x slabs, using distributed
filters and neighbor primitive gradients. Two/four-rank unequal partitions match
serial coefficients exactly; invalid local configuration/EOS data reject
collectively. Distributed dynamic fluxes, stepping and application dispatch
remain open. See [MPI coefficient-field evidence](validation/les-mpi-dynamic-field-20260918.md).

`compute_mpi_dynamic_les_rhs_3d` now uses distributed coefficients and two-layer
neighbor halos with the shared dynamic flux kernel. Owner-local RHS and the
global maximum diffusion estimate match serial results exactly in the focused
unequal-slab tests. Publication is collective and transactional. MPI dynamic
stepping and coupled/application dispatch remain to be connected; see
[MPI dynamic flux evidence](validation/les-mpi-dynamic-flux-20260918.md).

`sample_mpi_les_x_planes_3d` provides owner-routed raw donor planes for wide
regular MPI physical halos, with unequal/zero request counts and transactional
outputs ([communication checks](validation/mpi-les-donor-planes-20260925.md)).
`fill_mpi_les_boundary_ghosts_3d` now supplies the shared physical mapping,
repeated moving-wall transforms and pressure-target correction, including
one-cell local slabs ([halo checks](validation/mpi-les-physical-halos-20260925.md)).
`compute_mpi_dynamic_les_boundary_rhs_3d` and
`advance_mpi_dynamic_les_ssprk2_3d` now connect physical fluxes, global boundary
rates and SSPRK2 budgets ([operator checks](validation/mpi-dynamic-physical-20260925.md)).
Regular MPI spray application dispatch now supplies these physical inputs and
adds both LES exchanges to its boundary ledger
([wall/restart checks](validation/mpi-dynamic-spray-walls-20260925.md)).

`advance_mpi_dynamic_les_ssprk2_3d` now integrates owner-local periodic or physical-boundary dynamic
LES with collective dt agreement, global diffusion checks and transactional EOS
recovery after both stages. Two/four-rank variable-state checks match serial
state and temperature exactly. See [distributed stepping evidence](validation/les-mpi-dynamic-step-20260918.md).

MPI mixed-spray coupling and both MPI applications now accept the dynamic model
and filter controls. Controls agree collectively and enter restart identity;
periodic SGS half steps use distributed coefficient, flux and SSPRK2 APIs.
Same-rank exact continuation and 2-to-4-rank continuation pass for fixed/selected
H2O2. Regular MPI physical-boundary Dynamic LES now reaches the spray application,
with focused fixed-mechanism wall/restart evidence; the MPI AMR path
now supports physical boundaries and recursive stepping as described below. See
[MPI dynamic application evidence](validation/les-mpi-dynamic-application-20260918.md).

### Approximate LES filters

`les_test_filter` also accepts `box3`, `box5`, `gaussian3`, `gaussian5` and their
`opt` suffix variants. These use pinned polynomial formulas or optimized tables
for ratios 1–10, followed by the reference fallbacks. They retain signed weights;
filtered state may fail EOS validation, and very wide polynomial kernels may
fail the existing normalization guard. Serial/MPI use the same filter identities.
See [formulas, limits and focused checks](validation/les-filter-approximations-20260918.md).

Dynamic serial/MPI RHS and SSPRK2 APIs accept optional `harmonic_mean=.true.`
for coefficient/stress/thermal-quantity face interpolation. Same-sign values use
a signed harmonic mean; opposite signs and zero inputs produce zero. Omission
retains arithmetic interpolation. MPI rejects different interpolation settings
across ranks. Serial/MPI coupled APIs and fixed/selected applications expose
`les_harmonic_mean=.true.` with `les_model='dynamic'`; the default is `.false.`.
The setting enters checkpoint identity only when enabled, preserving the
arithmetic default context. Restart rejects a changed or removed true setting.
See [reference and harmonic parity checks](validation/les-harmonic-20260918.md)
and [application restart checks](validation/les-harmonic-application-20260918.md).

### Static LES isotropic stress

`smagorinsky_stress`, `wale_stress`, `vreman_stress`,
`compute_les_rhs_3d`, `compute_les_ghosted_rhs_3d` and
`advance_les_ssprk2_3d` accept optional trailing `ci` in `[0,1]`.
Omission or zero retains the deviatoric model. Positive `ci` adds the pinned
isotropic stress, including pressure work in the energy flux; it does not change
eddy viscosity or species/heat diffusion. It can act with `cs=0`.
Both SSPRK2 stages include its additional timestep estimate and EOS checks.
`advance_mpi_spray_les` also accepts trailing `ci`, checks its finite range and
rank agreement before stepping, and forwards it to both periodic and physical
boundary paths. Serial/MPI coupled APIs and fixed/selected application namelists
accept `les_ci` in `[0,1]`, default zero. Positive CI is allowed for Smagorinsky,
WALE and Vreman; dynamic LES derives its own isotropic coefficient and rejects
nonzero fixed CI. MPI coupling checks CI agreement and enables the static SGS
half steps with positive CI even when `les_cs=0`. Nonzero CI enters restart
identity; changing or removing it rejects restart. See the
[application checks](validation/les-isotropic-application-20260918.md).
See [formula and focused evidence](validation/les-isotropic-20260918.md) and
[MPI stepping checks](validation/les-mpi-isotropic-20260918.md).

### LES pressure outlets

The serial `compute_les_rhs_3d` and `advance_les_ssprk2_3d` library APIs accept
trailing optional `outlet_pressure(6)`, in Pa. Zero selects ordinary outflow;
positive pressure is allowed only on `outflow` faces. The shared characteristic
pressure-ghost routine supplies both halo layers before LES face gradients at
each stage. Smagorinsky, WALE and Vreman use the same boundary treatment,
including CI. The supplied-halo RHS continues to consume caller-prepared ghosts.
`advance_mpi_spray_les` now accepts the same outlet control. The shared
`fill_mpi_spray_ghosts` accepts optional `species` and `outlet_pressure`; positive
pressure requires species metadata. Pressure settings agree collectively, only
physical x-face owners apply those corrections, and tangential pressure ghosts
travel with the neighbor halos. Serial `advance_spray_reactive_3d` and
fixed/selected serial spray namelists now accept `outlet_pressure(6)`, passing
it to hydro, molecular transport and static LES. Positive pressure requires
outflow faces; nonzero settings enter checkpoint identity. Changed/removed
settings reject restart. MPI hydro and molecular-transport wrappers, mixed-spray
coupling, and fixed/selected MPI application inputs also accept the same control.
Pressure agreement is checked collectively before stepping. The shared physical
ghost routine supplies each operator's halo; nonzero settings enter MPI restart
identity. See [MPI application checks](validation/spray-mpi-pressure-outlet-20260918.md).
See [serial application evidence](validation/spray-pressure-outlet-20260918.md).
See [MPI pressure parity](validation/les-mpi-pressure-20260918.md).
See [focused evidence](validation/les-pressure-outlet-20260918.md).

### Static LES face fluxes for mesh coupling

`compute_les_ghosted_rhs_3d` optionally returns `fx`, `fy`, `fz` together.
For `cells=[nx,ny,nz]` and `nv` conserved components, use shapes
`(nv,0:nx,ny,nz)`, `(nv,nx,0:ny,nz)`, `(nv,nx,ny,0:nz)` respectively.
These are instantaneous flux densities in the positive coordinate directions,
without face-area or time factors. They are the same fluxes used for the RHS.
The caller supplies two halo layers, including any coarse/peer/physical donors.
All three outputs publish only on success; partial output sets or wrong shapes
are rejected. Boundary masks affect the ledger, not the returned flux values.
This supplies the fluxes used by the serial AMR LES integrator below.
See [focused evidence](validation/les-face-fluxes-20260918.md).

### Serial AMR static LES

The spray tree runtime exports `les_amr_patch_tree_mod` and
`advance_amr_patch_tree_les_3d`. It advances regular patch trees with common-step
SSPRK2, shared peer/coarse ghost filling, optional limited coarse donors,
physical boundaries including pressure outlets, conservative face synchronization
and reflux, averaging down and EOS recovery. It supports Smagorinsky, WALE and
Vreman with optional CI. Every patch/stage checks the existing LES diffusion
interval bound. State, temperature and boundary/correction diagnostics publish
together on success. `compute_amr_patch_tree_les_timestep_3d` selects a positive
common timestep from the same initial averaged state and ghost-aware diffusion
bounds. Supply a diffusion CFL fraction in `(0,0.2]`; its output is capped at
`1e6` seconds, including when all LES diffusion vanishes. It also returns the
largest patch diffusivity, which alone is not the grid-dependent interval bound.
Outputs remain unchanged on failure. The stepper still checks both stages and
can reject a selected interval after nonlinear evolution or EOS recovery.
The MPI spray tree runtime also exports
`les_mpi_amr_tree_mod::advance_mpi_amr_tree_les_3d`, using owner-local fields,
collective control agreement, distributed ghost queries, restriction and reflux.
It preserves all caller fields and diagnostics on collective failure. Empty
owners are supported. `compute_mpi_amr_tree_les_timestep_3d` now returns the
collective minimum patch interval and maximum diffusivity, with the same CFL
range/cap as the serial selector. It checks controls collectively, evaluates
owner-local ghosts and preserves both outputs on failure. LES subcycling
remains open. The MPI AMR spray application connects the common-step API below. See
[MPI timestep checks](validation/les-mpi-amr-timestep-20260918.md). See [distributed AMR evidence](validation/les-mpi-amr-tree-20260918.md).
See [timestep evidence](validation/les-amr-timestep-20260918.md).
EB grids are outside both APIs. The pinned PeleC reference also rejects
non-regular EB cells in all four LES closures; cut-cell LES would require a
separate model extension ([source audit](validation/les-reference-scope-20260925.md)).
See [AMR evidence](validation/les-amr-tree-20260918.md).

The serial stepper and timestep selector accept `model='dynamic'` on periodic,
prescribed reservoir inflow, extrapolated outflow and slip/adiabatic/isothermal wall trees.
Optional `test_model`, `test_ratio`, `grid_ratio`, `harmonic_mean` default
to `box`, 2, 1 and false. These controls are rejected for static models. Dynamic
coefficients replace the static `cs/prt/sct/ci` parameters.
Pressure-target outlets use the existing acoustic pressure extrapolation across
all donor layers. A wide halo can cross the physical domain boundary even for
an interior patch; correction uses that boundary rather than the patch edge.
Large pressure differences can make an outer donor inadmissible and reject the
transaction. See [pressure integration checks](validation/les-dynamic-amr-pressure-20260925.md).
Each stage builds sufficiently
wide donor halos, evaluates dynamic faces, synchronizes/refluxes and averages
down through the existing transaction. The diffusion bound uses the faces
contributing to each patch interior, including the boundary faces; unused outer
padding does not restrict the timestep.
The shared donor sampler now supports repeated moving-wall reflections for
halos wider than the domain, used by Dynamic wall
dispatch ([boundary-sampling evidence](validation/amr-wide-wall-ghosts-20260925.md)).
The ghosted Dynamic kernel also accepts local adiabatic/isothermal wall faces
and applies their heat flux before face-divergence publication
([kernel checks](validation/les-dynamic-wall-flux-20260925.md)). AMR selectors,
SSPRK2 and the coupled application forward local wall faces and temperatures;
interior/coarse-fine faces retain their donor fluxes. Moving walls retain the
existing molecular-viscosity requirement. See
[wall integration checks](validation/les-dynamic-amr-walls-20260925.md).
Reacting-spray/application dispatch is described below. See
[serial Dynamic AMR checks](validation/les-dynamic-amr-20260925.md).

The MPI selector and stepper have the same dynamic boundary dispatch and
filter controls. They agree on filter/model options before distributed donor
queries and retain owner-local fields and collective publication. Focused
two-rank serial/MPI state, temperature, timestep and ledger comparisons now pass,
including empty ranks and distributed ownership. See
[runtime record](validation/mpi-dynamic-amr-runtime-20260925.md).
Inflow/outflow boundary ledgers use the existing conservative reflux transaction;
see [open-boundary checks](validation/les-dynamic-open-20260925.md).

### MPI AMR reacting-spray LES coupling

The regular-grid `fill_reactive_boundary_ghosts_3d` helper now accepts 1..128
layers with repeated moving-wall reflection; default width remains two
([donor checks](validation/regular-wide-boundary-donors-20260925.md)).
`compute_dynamic_les_boundary_rhs_3d` and `advance_dynamic_les_ssprk2_3d` now
connect these donors to regular-grid Dynamic physical fluxes and boundary budgets
([operator checks](validation/regular-dynamic-boundaries-20260925.md)).
Regular serial spray application dispatch now includes these boundary exchanges
([application checks](validation/regular-dynamic-spray-20260925.md)). Regular MPI
spray also connects physical boundaries
([MPI application checks](validation/mpi-dynamic-spray-walls-20260925.md)).

The serial `advance_amr_patch_tree_les_level_3d` API advances one depth using
temporal coarse donors and returns SSPRK2 time-integrated faces for later
reflux. It preserves inactive levels and publishes outputs transactionally.
The serial `advance_amr_patch_tree_subcycle_les_3d` driver now schedules
ratio-squared child intervals, interpolates coarse endpoints and synchronizes
time-integrated fluxes after each child cycle. It returns per-level step counts,
composite boundary change and maximum correction/diffusivity transactionally.
`compute_amr_patch_tree_subcycle_les_timestep_3d` selects the coarse interval by
scaling each patch's diffusion rate with the cumulative squared refinement ratio
([selector checks](validation/les-subcycle-timestep-20260925.md)). Stage bounds
remain enforced during stepping.
The MPI timestep API accepts `subcycle=.true.` for the same coarse-interval
scaling, with collective agreement and owner-local diffusion evaluation
([MPI selector checks](validation/les-mpi-subcycle-timestep-20260925.md)).
`advance_mpi_amr_tree_les_level_3d` now supplies owner-local level stepping,
temporal coarse donors and integrated faces with collective rollback
([MPI level checks](validation/les-mpi-amr-level-20260925.md)).
`advance_mpi_amr_tree_subcycle_les_3d` connects ratio-squared recursive intervals,
owner-local coarse endpoints and accumulated flux synchronization. It publishes
state, temperature, boundary budget, maximum correction/diffusivity and optional
level counts only after the complete tree step succeeds
([MPI recursion checks](validation/les-mpi-amr-subcycle-20260925.md)).
The MPI spray application now dispatches this recursive driver when
`subcycling=.true.` ([coupling checks](validation/les-mpi-spray-subcycle-20260925.md)).
See [level integration checks](validation/les-amr-level-20260925.md) and
[recursive integration checks](validation/les-amr-subcycle-20260925.md).

`advance_mpi_amr_spray_reactive` accepts optional `les_settings` containing
`[eddy coefficient, PrT, ScT, CI]` and `les_model` (`smagorinsky`, `wale` or
`vreman`). Defaults disable LES. Two LES half steps surround the gas
reaction/transport/hydro split, inside the particle half steps. `subcycling=.true.`
selects recursive gas and LES stepping; false retains common stepping. The
application also selects the scaled LES diffusion bound in recursive mode.
Subcycling, coefficients and model must agree across ranks. Gas, particles, optional films, IDs and boundary/correction
diagnostics publish together after the entire coupled interval succeeds.
The fixed/selected MPI AMR spray application exposes these controls through
`&amr_spray_les`, binds nondefault settings to restart identity and applies the
initial LES diffusion estimate to its automatic timestep. Existing retries
handle later nonlinear rejection. See the
[input guide](../cases/mpi_amr_spray_3d/README.md#static-les) and
[application checks](validation/les-amr-application-20260918.md). See
[coupled evidence](validation/les-amr-coupled-20260918.md).

The coupled API also accepts `les_model='dynamic'` with periodic/inflow/outflow or slip/adiabatic/isothermal wall faces and optional
`les_test_model`, `les_test_ratio`, `les_grid_ratio`, `les_harmonic_mean`.
Dynamic mode runs with zero static coefficients. Effective settings are agreed
before particle work, and all state still publishes only on full success.
The application selects it in `amr_spray_les` and reads filter settings from a
separate `amr_spray_dynamic` group, preserving older static context serialization.
Fixed-mechanism continuous/restarted MPI application checks now pass, including
2-to-4-rank continuation and rejection of changed filter settings. Physical
validation remains outside this evidence
([runtime record](validation/mpi-dynamic-amr-runtime-20260925.md)).
Adaptive grid changes/rebalancing and selected H2O2 continuation also pass
[focused checks](validation/les-dynamic-adaptive-20260925.md); this does not
qualify arbitrary selected mechanisms or detailed-fuel ignition.

### MPI AMR uniform acceleration

`advance_mpi_amr_spray_reactive` also accepts `acceleration(3)` and
`source_change(nvar)`. A nonzero acceleration requires the source ledger;
components must be finite, bounded by `1e6` in magnitude and identical across
ranks. Two half kicks update gas and parcel velocities around the gas/LES
operators, preserving internal energy and recording external impulse/work.
Covered coarse cells receive the update but do not contribute twice to the
composite ledger. Stationary wall films remain attached to their walls.
Source diagnostics publish only with the entire accepted coupled interval and
remain separate from boundary exchange. The fixed/selected MPI AMR spray
application exposes `&amr_spray_force acceleration=... /` and retains cumulative
source inventory through changed-rank checkpoint/restart. See the
[input guide](../cases/mpi_amr_spray_3d/README.md#uniform-acceleration),
[application evidence](validation/spray-amr-force-application-20260918.md) and
[focused evidence](validation/spray-amr-force-20260918.md).

### MPI child ID ordering

The shared MPI computational splitting, KH shedding and wall-event splash APIs
assign new IDs by global parent ID and deterministic sibling order within each
collective event. Assignment uses compact parent/count metadata and publishes
with the particle transaction. Repartitioning the same event's parents therefore
does not change these IDs. Initial injection still uses its separate reservation
order; different physical event grouping can also change IDs. Future child IDs
can differ from older builds, while existing saved IDs remain unchanged on load.
See [computational splitting](validation/spray-split-id-order-20260918.md) and
[KH/splash checks](validation/spray-physical-child-id-order-20260918.md).
The level-set EB impact wrapper uses the same wall-event routine and has
[reversed-owner ID checks](validation/spray-eb-child-id-order-20260918.md).

### Experimental common-step SDC particle API

`advance_mpi_amr_mixture_coupling` accepts `motion_mode='kick_drift'`
(half-interval physical update, full-interval drift) or `'kick'`
(half-interval physical update, fixed position). The default `'coupled'` retains
its full-interval update. Exported `gas_change` remains an integrated source;
divide by half the supplied interval for a kick stage's rate.

`advance_mpi_eb_spray_sdc` connects these stages to common-step MPI EB gas SDC.
Each new-source evaluation starts from the same saved old-stage endpoint.
Only the final accepted kick is published, together with gas, temperature,
reaction history, particle IDs and inventories. Internal particle exchange is
separate from external gas forcing. This library path supports specular contacts;
specified-state/Mach inflow and isothermal-wall gas data reach the SDC operators.
Computational multiplicity splitting runs once after the accepted final kick,
with collective rollback on capacity failure and globally assigned child IDs.
Optional `particle_acceleration(3)` adds uniform external particle acceleration
(default zero), independently of the required gas force density `force` (N/m³).
Optional `gas_acceleration(3)` separately supplies gas acceleration (m/s²),
with momentum source `force + density*gas_acceleration` and local velocity work.
The recursive MPI EB gas driver also forwards this acceleration to every level
interval ([checks](validation/sdc-recursive-acceleration-20260927.md)); recursive
particle and film scheduling still requires integration.
The MPI EB spray application passes spray `acceleration` to both optional
acceleration arguments
([density-dependent forcing checks](validation/sdc-gas-acceleration-20260927.md)). Particle half kicks
bracket the coupled interval; each new kick is private to its SDC evaluation.
The final source inventory sums gas forcing and the accepted particle impulse
and kinetic work; `gas_particle_change` still contains internal exchange only.
The MPI EB spray application exposes common-step SDC and coupled restart
([application checks](validation/spray-sdc-app-20260927.md)). Stationary-film
thermal evolution is connected at library level below; film impact/application
coupling and recursive particle stepping remain open. See the [acceleration checks](validation/spray-sdc-force-20260927.md). See the
[focused evidence](validation/spray-sdc-kick-drift-20260927.md).

The boundary/splitting extension is recorded in the
[focused checks](validation/spray-sdc-boundaries-split-20260927.md).

### Film source export for SDC

`advance_mpi_amr_film_coupling` now accepts optional `gas_change`, an owner-local
field of time-integrated conserved gas feedback per fluid volume. It accumulates
accepted deposition contributions directly and publishes only with the complete
film/gas interval. Rejected trials and wall heat do not enter this field; wall
heat retains its separate ledger. Empty or disabled film evolution returns a
zero field. Presence must agree across ranks. The common-step SDC wrapper now
accepts `films`, `film_mass` and `film_settings` together, with optional
`film_steps`. Existing stationary films receive a half-interval thermal update
after each particle stage; each new evaluation restarts from the saved midpoint.
Only the final film state and old/final wall heat publish. Optional
`film_impacts=.true.` connects the existing EB deposition/splash model during
the old kick/drift stage and requires all film arguments. Its value must agree
across ranks. New stages do not move particles and cannot repeat the impact.
`model_rejected` reports empirical-model rejection; state and IDs still roll back.
The default retains specular particle contacts. The MPI EB spray application
now connects `wall_impacts=.true.` in common-step SDC, including film checkpoint
continuation ([application checks](validation/spray-sdc-impact-app-20260927.md)). See [impact checks](validation/spray-sdc-impact-20260927.md). See [SDC film checks](validation/spray-sdc-film-20260927.md).
See [focused evidence](validation/film-sdc-source-export-20260927.md).

### Recursive gas source interface

Serial and MPI EB recursive SDC drivers now accept paired `auxiliary_source`
and `auxiliary_change` arguments. The callback receives each level interval's
start time relative to the outer step, in addition to interval length, stage
and depth. The selected-level gas driver masks inactive levels and integrates
only physical leaves; the recursive driver accumulates accepted interval
inventories separately from external forcing. Callback-owned particle/film
state must remain private until the outer call succeeds. This is the gas-side
interface; recursive particle/film scheduling remains unimplemented. See
[focused checks](validation/sdc-recursive-auxiliary-source-20260927.md).

Recursive gas drivers also accept `route_auxiliary=.true.` to retain source
feedback on levels other than the one being advanced. Only uncovered fluid
cells on those levels contribute. Old/final-new increments wait in private
buffers until the interval at `min(source depth,target depth)` closes; they
are then applied after flux synchronization, followed by EOS recovery and
restriction. The auxiliary inventory includes these increments once.
The default keeps the previous selected-level masking. This supplies source
routing, not the remaining per-particle recursive schedule. See
[cross-level checks](validation/sdc-cross-level-source-20260927.md).
