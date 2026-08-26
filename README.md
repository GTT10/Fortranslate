# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The `0.152.0` milestone contains the serial verification suite, eight optional
MPI executables, and runnable serial and sparse-MPI one-dimensional
reactive AMR applications with solution-driven dynamic regridding and
molecular transport. The sparse MPI driver can write an intermediate
patch-tree checkpoint and restart it with a different MPI rank count. The
serial two-dimensional EB AMR driver can create, move, resize, remove, and
re-create one fine rectangle from temperature-gradient tags while preserving
its composite conserved state, and can compose active-cell chemistry with the
two-level or root-only EB hydrodynamic path. Formatted checkpoints preserve
that lifecycle state for a later serial restart. With multipatch mode enabled,
the same application clusters disconnected tags into multiple fine rectangles,
selects a set-wide CFL limit, advances their reactive hydrodynamics and
chemistry transactionally with one flux register per child, periodically
rebuilds the set, and writes one CSV per active child. A dedicated formatted
checkpoint preserves the complete ordered patch set for transactional serial
restart. A configured static fine rectangle may also meet an outflow physical
boundary: its physical-side exterior state is extrapolated from the current
fine boundary cell while the remaining coarse/fine sides retain coarse-time
interpolation and conservative reflux. Temperature-gradient tagging and both
single- and multipatch planners can now create such outflow-side rectangles
dynamically, including one-sided gradient detection on root boundary cells.
The single-patch, separated sibling-patch, and strictly nested three-level
lifecycles also accept mixture molecular transport. Every parent/child pair
advances SSPRK2 transport with ratio subcycling, time-interpolated parent
exterior states, an independent open-area diffusive flux register, reflux,
and hierarchy-wide average-down inside the transactional `R-T-H-T-R`
composition.
The MPI EB AMR bridge executes chemistry, owner-tiled finite-band root hydro
and transport, and fine ratio subcycles on deterministic physics owners. Fine
transport owners assemble patch-plus-two start, uncorrected-end, and current
corrected state/temperature directly from intersecting root tile owners, along
with child-intersecting x/y flux fragments. They reconstruct the four-edge
exterior context, assemble the compact coarse register, perform reflux locally,
and return corrected support directly to those tile owners. Fine state remains
on its owner, and no complete root state or flux array enters the child route.
Collective input consensus and deferred publication preserve exact all-rank
rollback.
One outer MPI transaction now composes those owner operators as
`R-T-H-T-R`. It publishes the root, children, limiter minimum, and chemistry,
hydro, and transport counters only after every stage succeeds; a rejection in
the middle hydro stage discards valid reaction and transport prefixes.
The 2D EB bridge also has a rank-local sparse payload container. Each rank
allocates conserved state and temperature only for its owned root tiles and
children; an explicit materialization boundary reconstructs the replicated
hierarchy when a legacy operator or output path still requires it. A second
materialization boundary gathers each root tile and child only to a
caller-selected root rank. Non-root ranks keep the complete output unallocated,
and one packed point-to-point message is sent per remote entity. This is the
field-sparse foundation for checkpoint and output adapters. Sparse MPI wrappers
now connect that boundary to the existing formatted multipatch checkpoint and
root/child CSV writers. Only the selected writer allocates the complete fields
or touches the files, and its I/O result is returned collectively. The inverse
restart boundary reads the formatted checkpoint only on a selected root and
sends each root tile or child directly to its current owner. Non-root complete
read buffers stay unallocated and no numerical-field broadcast is used. The
restart caller now supplies a replicated geometry-only child descriptor; it no
longer carries replicated child state or temperature fields. Compatibility
wrappers retain the former full patch-set API for existing callers.
Separately, a geometry-only 2D EB patch-tree topology now represents any number
of refinement relations. Each level may branch across multiple ordered parents
and children, and a complete replacement plan commits only after parent links,
EB geometry, patch nesting, and sibling separation validate.
Chemistry now runs directly on those sparse owner allocations. Root tiles and
children are reacted locally with covered cells masked; only post-reaction
average-down communicates numerical state. Each child owner sends one
coarse-footprint restriction buffer only to the distinct root tile owners that
intersect it; unrelated ranks receive nothing. Those owners recover temperature
and apply the restriction locally. Chemistry no longer materializes a complete
hierarchy or broadcasts child restrictions.
The complete sparse `R-T-H-T-R` transaction now remains on sparse owners from
input through commit. Both chemistry half-steps, both SSPRK2 transport
half-steps, and hydro call their direct sparse entrypoints; no complete child
hierarchy is materialized between operators. Counts and the limiter minimum are
still published only with the final outer commit.
Hydro also has a direct sparse entrypoint. Each root tile owner receives only
the neighboring row fragments required for a six-row halo, advances its own
bounded EB band, and retains its stage-start, stage-end, current corrected
state/temperature, and uniquely owned flux rows locally. Each child owner
assembles patch-plus-two start/end/corrected state directly from intersecting
tile owners and extracts the four-edge interpolation context locally. The same
tile owners route their x/y interface-flux fragments directly to the child.
After compact register accumulation, ratio subcycling, and reflux, corrected
support returns directly to those tile owners in deterministic child order.
The final corrected root rows commit locally. No complete hydro state,
temperature, or flux result is assembled on a root physics owner, and there is
no post-compute tile-result or final-scatter traffic.
SSPRK2 molecular transport now follows the same sparse ownership boundary.
For each Euler stage, a root tile owner receives only the neighboring row
fragments needed for its six-row transport/StateRedist guard, advances that
target EB band, and retains its owned stage-start, stage-end, corrected state,
temperature, and unique flux rows locally. A periodic y-boundary tile uses a
boundary-anchored cyclic band built from two contiguous source-row fragments.
One extra row isolates
the required six-row guard from the deliberate internal gap; a small root that
cannot hold that guard uses the complete root band as its local compute band.
No post-compute complete root state, temperature, or flux bundle is assembled.
Each fine child assembles patch-plus-two stage-start,
uncorrected-end, and current corrected state/temperature directly from
intersecting root tile owners. Those owners also send their coarse x/y flux
fragments. The child extracts its exterior context, assembles its register,
performs ratio subcycling, fine-flux accumulation, and reflux locally, then
returns corrected fragments directly to the same root tile owners. The final
SSPRK2 root blend and EOS
temperature recovery are cell-local on those owners, so that step needs no root
gather or scatter. For cut-interface closure, tile owners sum their physical-
boundary flux contributions and combine one conserved `nvar` vector with
`MPI_Allreduce`; unrefined-cell corrections remain tile-local. A late child
failure leaves every sparse allocation bitwise unchanged and publishes zero
work or transfers.
EB flux registers now store only the patch-plus-one-cell correction support.
For sparse transport, each root tile owner retains its computed x-flux rows and
uniquely owned y-faces, then sends only child-intersecting fragments directly
to that child owner. The child assembles patch-local face rectangles and
accumulates its coarse register locally; the root physics owner no longer
builds or sends that register. The established complete-root accumulation call
is a compatibility wrapper over the same compact kernel. Remote child owners
allocate no complete root state, temperature, or flux field; fine state and
coarse correction no longer make a root-physics-owner round trip for reflux.
Exterior-context extraction also accepts globally indexed patch-plus-one
coarse start/end state and temperature support. The complete-root entrypoint
is a wrapper over that support kernel. Sparse MPI uses it on the child owner
after direct tile-fragment assembly without changing interpolation or EOS
recovery.
The sparse hierarchy also selects its own stable coarse interval. Every root
tile evaluates its EB hydro and molecular-transport limits directly on its
owner, while fine children do the same and scale their stable fine steps by the
refinement ratio. One communicator-minimum reduction selects the interval with
zero root-field traffic. Invalid owner state rejects collectively with zero
published dt and transfer count.
The same sparse hierarchy now owns a public full-physics time loop. It
recomputes that distributed stable interval before every `R-T-H-T-R` step,
clips the last interval exactly to the requested target time, and publishes
time, step, operator-count, limiter, and timestep-traffic diagnostics only
after each complete step commits. A later step limit or physics rejection
preserves the already committed prefix and its matching diagnostics.
An explicit sparse regrid can now change the ordered child-patch topology,
reapply the serial EB overlap/prolongation rules, rebuild deterministic
subcycle-weighted ownership, and return to one-copy sparse storage. The whole
  distribution/template/state replacement is transactional. Regrid now remains
  field-sparse: old child owners restrict directly to intersecting root owners,
  root tiles initialize each distinct new child owner, and retained overlap
  moves directly from its old owner to its new owner.
The sparse public clock can also evaluate temperature-gradient tags at a
caller-selected accepted-step cadence. Only the root physics owner constructs
the patch plan; compact topology metadata is broadcast, while the caller's
geometry callback rebuilds each planned EB child. Physics and any due regrid
commit as one transaction, so a geometry or topology failure leaves the step,
clock, hierarchy, counters, and published traffic unchanged.
The EB transfer foundation also accepts one strictly nested three-level
hierarchy, computes its composite integral without double counting, and
average-downs its generic or reactive state from the deepest level to the root
as one rollback-safe transaction.
The same hierarchy can advance reactive EB hydrodynamics recursively with one
root update, ratio-subcycled middle and finest updates, an independent flux
register at each interface, and final deepest-first synchronization. When the
finest interface crosses the embedded boundary, a multilevel conservation
closure returns the measured mass, total-energy, and species residual to
uncovered active middle cells before the outer reflux.
The three-level driver now composes active-cell reaction half-steps on the
root, middle, and finest meshes around that recursive hydro transaction. A
second deepest-first reactive average-down makes the post-chemistry hierarchy
authoritative before all three levels are published together.
With `three_level_enabled`, the public EB AMR application constructs that
static hierarchy from two nested namelist rectangles, selects a root timestep
from all three active-cell CFL limits, advances to the requested final time,
and writes separate root, middle, and finest CSV fields. A dedicated
three-level checkpoint preserves all three conserved and temperature fields,
accepted time and step accounting, and the full nested topology for a
transactional serial restart. With dynamic regridding enabled, the middle
level remains fixed while EB-aware temperature-gradient tags can move and
resize the finest patch transactionally after initialization and accepted
steps. A distinct dynamic three-level checkpoint stores that committed finest
rectangle, regrid count and compatibility controls so a split run resumes the
actual hierarchy and cadence rather than the configured seed.
The single-level reactive EB application now also composes mixture viscosity,
thermal conduction, species diffusion, and barodiffusion with chemistry and
hydrodynamics. Diffusive Cartesian-face fluxes are interpolated to EB face
centroids, weighted by open area, stabilized with EB-aware species limiting
and StateRedist, and use an adiabatic slip/impermeable embedded-wall closure.

### `pelef`: one-dimensional Euler solver

- constant-`gamma` ideal-gas EOS;
- PCM, componentwise primitive PLM, and time-traced characteristic PLM;
- order-2 or PeleC-style five-point order-4 limited slopes;
- optional pressure/velocity shock flattening;
- Rusanov or qualified single-species PeleC-style Riemann fluxes;
- outflow and periodic boundaries;
- Sod, Shu-Osher, and symmetric planar Sedov-type regressions.

### `pelef2d`: two-dimensional Euler scaffold

- uniform periodic Cartesian mesh;
- x/y directional Riemann fluxes through explicit momentum rotation;
- limited primitive slopes and normal characteristic tracing;
- CTU-style transverse half-step corrections with positivity scaling;
- one unsplit conservative update;
- periodic isentropic-vortex analytical and convergence regressions.

### `pelef_ms`: passive multispecies Euler solver

- runtime species count with conserved `rho*Y_k` components;
- checked mass-fraction conversion, positivity, and species closure;
- species fluxes tied exactly to the shared mass flux;
- one-dimensional characteristic tracing and two-dimensional CTU transport;
- MultiSpecSod and periodic species-wave regressions.

This older passive path intentionally retains the constant-`gamma` hydro baseline.

### `pelef0d`: thermodynamics and toy-reactor verification

- species molecular weights and NASA7 thermodynamic polynomials;
- mass-based mixture molecular weight, gas constant, `cp`, `cv`, `gamma`, enthalpy, internal energy, and frozen sound speed;
- ideal-gas pressure/density conversion;
- bracketed Newton/bisection inversion from specific internal energy to temperature;
- a synthetic constant-volume two-species isomerization reactor;
- isothermal analytical and adiabatic energy-conservation gates.

### `pelef0d_h2o2`: generated elementary H2/O2 kinetics

- runtime elementary-reaction records with arbitrary reactant/product stoichiometry;
- reversible Arrhenius rates and NASA7 equilibrium constants;
- molar concentrations, progress rates, production rates, and mass-fraction source terms;
- JSON-to-Fortran mechanism generation with a committed-source cleanliness gate;
- an adaptive explicit RK4 constant-volume, adiabatic reactor;
- a seven-species, four-reaction H2/O2/N2 subset selected from Cantera `h2o2.yaml`;
- live trajectory and exact-state production-rate comparison against Cantera 3.2.


### `pelef_transport_probe`: dilute-gas transport verification

- Lennard-Jones transport records for H2, H, O, O2, OH, H2O, and N2;
- Chapman--Enskog pure viscosities and binary diffusion coefficients;
- Wilke mixture viscosity;
- modified-Eucken pure conductivity and Mathur mixture conductivity;
- mixture-averaged species diffusion coefficients;
- live qualification against Cantera 3.2 at four temperature, pressure, and
  composition states.

This is a deliberately qualified ideal-gas subset. It does not yet reproduce
PelePhysics polynomial transport fits, polar corrections, Soret diffusion,
multicomponent diffusion, or bulk viscosity.

### `pelef_reactive_1d`: general-EOS reactive Euler solver

- conserved state `(rho, rho*u, rho*v, rho*w, rho*E, rho*Y_k)`;
- NASA7 composition-dependent pressure, temperature, heat capacities, ratio of specific heats, and frozen sound speed;
- safeguarded conserved-to-primitive recovery through `e(Y,T) -> T` inversion;
- PCM, frozen-composition characteristic PLM, monotone primitive PPM, or
  time-traced frozen-composition characteristic PPM reconstruction;
- selectable Rusanov, general-EOS HLLC, or NASA7 PeleC-style acoustic flux
  with exact species-flux closure;
- SSPRK3 time integration for the semidiscrete primitive-PPM path;
- PeleC-style parabolic profile integration over the `u-c`, `u`, and `u+c`
  waves for characteristic PPM;
- optional PeleC one-dimensional shock flattening and bounded
  Colella--Woodward contact steepening on characteristic PPM;
- periodic or outflow boundaries;
- cell-local adiabatic constant-volume chemistry;
- optional shear viscosity, Fourier heat conduction, mixture-averaged species
  diffusion, barodiffusion, correction velocity, and species-enthalpy flux;
- explicit SSPRK2 transport with a parabolic timestep gate;
- symmetric reaction--transport--hydro--transport--reaction splitting;
- homogeneous-reactor reduction, smooth density/composition-wave convergence,
  discontinuous material-contact resolution, primitive/characteristic PPM
  convergence, periodic strong-shock flattening, and nonuniform
  reactive-hotspot regressions.


### `pelef_reactive_2d`: general-EOS reactive CTU solver

- uniform periodic Cartesian mesh with the same conserved reactive state as the 1D path;
- composition-dependent NASA7 pressure, temperature, heat capacities, and frozen sound speed;
- directional general-EOS Rusanov, HLLC, or PeleC-style acoustic fluxes
  through explicit momentum rotation;
- PCM, frozen-composition characteristic PLM, or time-traced
  frozen-composition characteristic PPM in both coordinate directions;
- optional bounded contact steepening and PeleC-style shock flattening on the
  characteristic-PPM normal predictor;
- provisional face fluxes, conservative CTU transverse half-step corrections, and EOS-based positivity scaling;
- species, momentum, and total energy corrected together so `sum(rho*Y_k)=rho` remains coupled to the hydro update;
- cell-local chemistry with symmetric reaction--transport--hydro--transport--reaction splitting;
- optional x/y Newtonian viscosity, Fourier conduction, mixture-averaged species diffusion, barodiffusion, correction velocity, and species-enthalpy flux;
- exact diagonal density/composition-wave convergence, x/y one-dimensional
  reduction, material-contact sharpening, oblique strong-shock flattening,
  periodic vortex, and reacting-hotspot regressions.

The reactive 2D path supports matched periodic pairs, slip/no-slip walls,
adiabatic/isothermal wall temperatures, fixed-state inflow, and zero-gradient
outflow. Walls are species-impermeable by default or may impose a prescribed
zero-net-mass species-conversion flux with consistent species-enthalpy
transport. Periodic cases retain the qualified characteristic-PPM plus CTU
path; physical faces use boundary-aware ghost reconstruction, an exact
impermeable inviscid wall flux, and boundary-aware molecular transport.

The reactive applications can select either the verified seven-species,
four-reaction elementary subset or the full ten-species, 29-reaction H2/O2
mechanism with third-body, falloff, Troe, and adaptive implicit chemistry. The
characteristic projection is a qualified frozen-composition ideal-gas-mixture
approximation, not full PeleC/PelePhysics general-EOS characteristic parity.

### Embedded-boundary geometry foundation

`eb_geometry_2d_mod` converts a nodal level set into bounded Cartesian
cell-volume fractions, x/y face-area fractions, and regular/cut/covered cell
types. Positive level-set values define fluid. Each cell uses two affine
triangles, so planar interfaces are exact and curved interfaces converge under
refinement. Cut cells also carry the physical embedded-boundary length,
centroid, and a unit normal directed from solid toward fluid. The reactive EB
kernel recovers the general-EOS pressure and applies the
stationary impermeable slip-wall momentum flux on an arbitrarily oriented cut
wall. It also converts the integrated wall flux into a volume-normalized source
for each cut cell without mass, energy, or species leakage. Curved-cell force
uses the integrated interface-normal vector rather than multiplying total
length by an averaged unit normal. The conservative divergence combines this
wall contribution with shared Cartesian fluxes weighted by open-face fraction,
and preserves a uniform stationary pressure field for planar and circular
embedded boundaries. Small-cell time integration now has a conservative
first-order FluxRedist path: it blends a cut-cell update with its
volume-weighted face-connected neighborhood, redistributes the removed
extensive update, and commits a forward update only after every active reactive
state passes EOS recovery. Weighted StateRedist forms normal-directed
neighborhoods with target volume fraction `0.5`, accounts for cells shared by
overlapping neighborhoods, conserves every volume-weighted state component,
and applies the same transactional EOS gate. Its default `max_order=0` path
retains neighborhood averages. Selectable `max_order=2` stores normalized
fluid-volume centroids, fits limited linear `Qhat` slopes over active
neighborhoods, evaluates them at every merge recipient, and adds a
conservation-compatible recipient maximum-principle limiter.
A selectable PCM or frozen-composition characteristic-PLM EB hydro path now
constructs reactive Riemann fluxes only on open Cartesian faces, uses
zero-gradient domain faces, linearly interpolates face-center fluxes to the
open-face centroid, combines them with the integrated slip-wall pressure
force, and completes the step through weighted StateRedist and EOS recovery.
PLM slopes use only two-sided active-cell stencils; cells beside covered or
outer cells fall back locally to zero slope. The `pelef_reactive_eb_2d`
application now reads plane or circular geometry from a namelist, initializes a
general-EOS multispecies state, advances to a requested final time with an
active-cell CFL limit, reports volume-weighted diagnostics, and writes cell
geometry and primitive fields to CSV. This qualified runnable path has
zero-gradient outer faces and optional active-cell chemistry. Its
Strang sequence applies half reactions only to active cells and can insert
SSPRK2 molecular-transport half steps around the EB hydro transaction while
leaving covered cells bitwise unchanged. Transverse reconstruction settings
remain rejected instead of silently ignored.

The EB and AMR foundations now meet at one qualified static two-level transfer.
An aligned rectangular fine patch is restricted with fine fluid-volume weights,
the corresponding composite integral counts uncovered coarse cells and fine
cells exactly once, and reactive restriction recovers every active parent
temperature before committing either state or temperature. Covered parents
retain their original reactive data. This is a serial synchronization kernel,
not yet a time-advancing EB AMR application.

The same static hierarchy now owns an EB-aware flux register. Coarse and fine
steps accumulate time-integrated face fluxes independently; open-face fractions
and physical subface lengths produce one correction on each exterior coarse
cell. A regular cell receives that correction directly. A cut cell keeps its
fluid-volume share and redistributes the remainder over connected 3-by-3
neighbors, with any share landing below the fine patch transferred to all of
that parent's fine children. Reactive re-reflux commits both levels and both
temperature fields only after every active cell passes EOS recovery.

The static hierarchy can now initialize its fine rectangle by piecewise-
constant coarse injection and advance reactive EB hydrodynamics over one
coarse interval. The coarse level takes one step while the fine level takes
`r` steps. Each fine-patch boundary face obtains a time-interpolated conserved
state from the adjacent coarse cell, followed by EOS temperature recovery.
The fluxes that actually advanced both levels feed the EB flux register;
re-reflux and reactive average-down then synchronize the hierarchy in one
transaction.

`pelef_reactive_eb_amr_2d` makes that hierarchy runnable from one input file.
It builds a strictly internal fine rectangle, initializes it from the coarse
state, selects each coarse timestep from both level CFL limits, advances until
the clipped final time, and writes separate synchronized coarse and fine
geometry/state CSV files. Optional solution-driven regridding tags active
coarse cells by relative and absolute temperature jumps, buffers their bounding
rectangle, averages the old fine patch down, injects the new patch from coarse
data, and retains every overlapping same-resolution fine cell exactly.
When configured to remove an untagged patch, the driver conservatively collapses
the child into the root, releases its arrays and geometry, advances with the
single-level CFL and hydro path, and re-creates the patch by PCM when tags
return. Fine CSV output is omitted while the child is inactive.

With chemistry enabled, each active AMR interval applies a reaction half-step
on both levels, the existing subcycled EB hydro/reflux transaction, a second
reaction half-step, and final fine-to-coarse average-down. Covered cells are
masked from chemistry. A root-only lifecycle interval reuses the qualified
single-level EB Strang path, and any chemistry, hydro, or EOS failure leaves the
complete hierarchy unchanged.

The EB AMR input can request periodic checkpoints, stop immediately after a
successful write, or restart from a prior file. A checkpoint records the
coarse state, optional fine state and current patch bounds, time/step/regrid
metadata, minimum accepted timestep, base density, species ordering, and the
physics and regrid settings required for compatible continuation. Restart
rebuilds both EB geometries from the input, recovers active-cell temperatures
from the conserved state, and publishes the hierarchy only after the complete
file and end marker pass validation. Final time, maximum steps, output paths,
and checkpoint controls may change between runs.

The two-level EB kernel can also cluster disconnected tag components into a
deterministically ordered collection of separated fine rectangles. Buffer and
minimum-size expansion are applied per component, nearby candidates are
coalesced to preserve the 3-by-3 redistribution neighborhood contract, and
clamped to the complete root domain. Boundary cells use available one-sided
temperature differences, so a single- or multipatch plan may meet any outflow
side. A patch-set transaction
supports PCM creation, exact old/new fine-overlap retention, conservative
average-down, composite integration, movement, repartition, and removal.
Hydrodynamics advances the root once, subcycles every child with its own EB
flux register, then refluxes and averages down all children. The matching
Strang transaction applies active-cell chemistry on the root and every child
around that hydrodynamic interval and restores the complete hierarchy after
any rejected stage.

Setting `multipatch_enabled = .true.` connects that patch set to the public
`pelef_reactive_eb_amr_2d` lifecycle. The input exposes the maximum tag gap,
the initial and periodic planners can replace the configured seed rectangle
with zero or more children, the timestep is the minimum root/child CFL limit,
and accepted steps retain the existing regrid cadence and counters. Fine CSV
paths receive deterministic `_patch0001`, `_patch0002`, ... suffixes.

Multipatch mode has a separate versioned formatted checkpoint schema. It stores
the root, ordered child count, every child's actual bounds and state, accepted
time and timestep, step/regrid counters, base density, species order, and a
strict physics/topology compatibility signature. Restart rebuilds every EB
geometry, recovers active temperatures through the EOS, validates the complete
set and end marker in private candidates, and only then publishes the restored
hierarchy. The earlier single-patch schema and its inputs remain unchanged.

Setting `three_level_enabled = .true.` instead constructs one static middle
rectangle from root indices and one finest rectangle from middle indices. The
finest rectangle retains a two-cell middle margin. The public timestep is the
minimum root-equivalent stability limit from all three levels, and every
accepted interval uses recursive subcycling, independent interface registers,
the EB-cut conservation closure, active-cell Strang chemistry, and final
deepest-first synchronization. Successful completion writes distinct root,
middle, and finest CSV files. Scheduled and final three-level checkpoints use
a dedicated schema and may stop and resume the same hierarchy without changing
the established single-patch or patch-set formats. Three-level mode remains
mutually exclusive with multipatch siblings. Its dynamic path keeps
the middle patch fixed, retains the finest patch, ignores tags outside its
two-cell-safe planning region, and does not yet support checkpoint/restart.

Unsplit transverse prediction, fourth-order StateRedist slopes,
periodic/ghost-cell neighborhoods, thermal/catalytic wall physics,
coarse-to-fine spatial slopes, multipatch EB AMR molecular
transport, dynamic middle/root topology, arbitrary-depth EB levels,
transport-enabled checkpoint/restart, arbitrary-depth EB physics recursion,
and MPI distribution are not yet
connected. The public EB AMR application now owns
either restartable sibling rectangles or an
explicit three-level hierarchy with an optionally dynamic finest patch.

The separate EB patch-tree core now owns reactive conserved state and
temperature on arbitrary-depth, branching topology. It initializes children
from their actual parents, synchronizes deepest-first, evaluates the complete
composite conserved vector, and transactionally migrates fields through a
whole-tree topology replacement. Same-resolution physical overlap is retained
only after local EB geometry checks; EOS or conservation failure leaves the
accepted tree unchanged. This core is not yet connected to the public physics
time loop or MPI ownership.

The replicated MPI-owner EB AMR hydro path now decomposes the root update over
its distributed y-tiles. Each tile owner advances a bounded six-row halo band,
publishes only its owned cells and uniquely assigned faces, and participates in
a collective root assembly. The previous single root-owner full-level advance
and four full-root broadcasts are absent from this path. Fine children retain
their established owner subcycling and deterministic reflux order. The sparse
root path remains a later conversion boundary.

### MPI one-dimensional verification

With `PELEF_ENABLE_MPI=ON`, five domain-decomposed executables verify:

- uneven non-replicated block decomposition and periodic halo exchange;
- conservative multispecies Euler transport;
- distributed adaptive implicit full-H2/O2 chemistry scheduling;
- general-EOS molecular transport with viscosity, conduction, barodiffusion,
  correction velocity, and species-enthalpy transport;
- transactional reaction--transport--hydro--transport--reaction splitting;
- ordered gather output, global timestep/conservation reductions, and
  complete-field parity for 1, 2, 4, and 8 ranks.

A sixth MPI executable exercises sparse AMR distribution. Compact hierarchy
and owner metadata are replicated, but each root/fine field payload exists
only on its deterministic work-weighted owner. Patch work may use raw cell
count, hyperbolic `r` subcycling, or parabolic `r^2` subcycling; explicit and
tag-driven regrids preserve the selected model. Rank-local hyperbolic and
parabolic stability limits are reduced to one communicator-wide coarse-step
limit without gathering patch fields. Point-to-point transfers cover
same-level halos, parent/child ghost data, boundary fluxes, shared-face
corrections, average-down, explicit regrid prolongation, retained overlap, and
owner migration. Chemistry, recursive hydro, parabolic `r^2` transport, and
the transactional `R-T-H-T-R` interval run on owners alone. Both explicit and
solution-tagged topology rebuilds remain field-sparse; owner-local tagging
shares only compact integer plan metadata. The 1/2/4/8-rank gates compare the
gathered hierarchy, fields, ghosts, counters, conservation, and rollback with
the serial patch-tree implementation.

The seventh executable, `pelef_mpi_amr_reactive_1d`, is the public sparse AMR
driver. It reads the reactive namelist, builds owner-local tagged patch trees,
selects distributed hydro/transport timesteps, advances `R-T-H-T-R`, regrids
at the requested cadence, and writes an ordered composite AMR CSV. Persistent
field payloads remain globally single-copy during the time loop; a complete
tree is materialized for final diagnostics/output and scheduled checkpoints.
Checkpoints store no owner map: restart rebuilds deterministic ownership for
the active communicator, allowing a two-rank run to resume on four or eight
ranks.

The eighth executable, `pelef_mpi_eb_amr_patch_2d`, establishes the first
two-dimensional EB AMR distribution boundary. It partitions the root into
contiguous y-tiles, assigns separated fine siblings as independent entities,
and balances raw, hyperbolic, or parabolic subcycle-weighted work with 64-bit
accounting. Geometry and patch metadata remain replicated, while root-tile and
child state/temperature payloads are synchronized from one authoritative
owner. Its 1/2/4/8-rank gate checks collective topology agreement, exact
ownership accounting, authoritative payload recovery, and transactional
rejection of invalid or inconsistent work models. The same gate now advances
active-cell chemistry exactly once on each root-tile or child owner, reaches
serial patch-set parity after owner broadcasts and fine-to-root average-down,
and rolls every rank back after a late owner-side reactor rejection. Reactive
EB hydrodynamics now advances the complete root level once on its exclusive
physics owner and advances every child on its patch owner with ratio
subcycling, owner-local flux registers, sequential reflux, and final
fine-to-root average-down. Its result matches the serial multipatch path and a
late child-owner failure leaves every rank unchanged.

### One-dimensional AMR

The AMR layer provides:

- a static two-level hierarchy with an interior refined patch;
- an arbitrary-depth hierarchy foundation built from nested level interfaces;
- arbitrary component counts and integer refinement ratios;
- MC-limited, conservative piecewise-linear prolongation;
- volume-average restriction and covered-cell synchronization;
- refinement-ratio level subcycling;
- time-integrated coarse/fine flux registers and reflux;
- a composite-integral gate proving conservation across both patch interfaces;
- component-selectable normalized-gradient tagging with absolute noise floors;
- buffered, minimum-width single-patch planning with explicit boundary checks;
- conservative patch creation, movement, resizing, and removal;
- exact retention of same-resolution fine data where old and new patches overlap;
- reactive general-EOS states advanced on both coarse and fine levels;
- refinement-ratio fine-level subcycling with time-interpolated coarse ghosts;
- coarse/fine flux accumulation, reflux, and average-down every coarse step;
- symmetric chemistry--hydro--chemistry splitting across the AMR hierarchy;
- transactional hierarchy rollback and periodic solution-driven regridding;
- ordered composite CSV output with exact domain-coverage checks;
- optional limited primitive-variable PLM with SSPRK2 time advancement;
- time-averaged SSPRK2 interface fluxes used consistently for reflux;
- four-layer physical and coarse/fine ghost storage for PPM stencils;
- primitive or time-traced characteristic PPM with SSPRK3 level advancement;
- optional PeleC-style hybrid WENO5-JS, WENO5-Z, WENO7-Z, or WENO3-Z edge
  reconstruction inside the characteristic-PPM path;
- SSPRK3 effective interface fluxes used consistently for reflux;
- conservative, MC-limited parent-to-fine ghost interpolation at subcycle
  midpoint times;
- AMR viscosity, Fourier conduction, and mixture-averaged species diffusion;
- parabolic fine subcycling with time-interpolated coarse transport ghosts;
- diffusive flux-register reflux and covered-cell average-down;
- symmetric reaction--transport--hydro--transport--reaction composition;
- cumulative multilevel subcycle schedules and deepest-to-root synchronization;
- a four-level, mixed-ratio conservation gate across every interface;
- arbitrary-depth reactive state and temperature ownership;
- recursive hydro subcycling and diffusive `r^2` subcycling at every depth;
- recursive chemistry and reaction--transport--hydro splitting;
- a three-level reactive conservation, positivity, closure, and synchronization
  gate;
- a runtime `amr_max_levels` limit with tag-driven nested hierarchy creation;
- periodic conservative multilevel hierarchy rebuilds during simulation;
- recursively ordered, exact-coverage multilevel composite CSV output;
- a runnable three-level hotspot case and structural output gate;
- physical-coordinate overlap transfer across changed multilevel hierarchies;
- exact retention of aligned old fine state and temperature data;
- nested fine patches touching an outflow physical boundary;
- physical-side PPM/WENO ghost fill with reflux restricted to the remaining
  coarse/fine interface;
- configured two-level EB fine patches touching an outflow physical boundary,
  with zero-gradient physical-side exterior state and no physical-side reflux;
- ordered sets of disjoint fine patches over one parent level;
- disconnected-tag clustering with deterministic buffer/minimum-width
  expansion and automatic coalescing of adjacent candidates;
- set-wide conservative prolongation, average-down, reflux, and composite
  integration without double counting covered parent cells;
- conservative patch-set movement, repartition, removal, and exact retention
  of same-resolution fine overlap across old/new patch intersections;
- fixed two-level reactive WENO7-Z hydro subcycling on two separated fine
  patches with per-patch flux registers;
- fixed two-level multipatch reaction--transport--hydro splitting with `r^2`
  diffusive subcycling, per-patch reflux, and transactional rollback;
- tag-driven two-level multipatch creation, movement, repartition, removal,
  overlap retention, runtime statistics, and ordered composite CSV output;
- arbitrary-depth separated patch trees with explicit parent ownership,
  mixed per-level refinement ratios, recursive conservative prolongation,
  deepest-to-root average-down, and exact composite integration;
- parent-owned patch-tree flux registers with transactional deepest-to-root
  reflux and covered-cell synchronization across every branch;
- static arbitrary-depth reactive patch trees with recursive PCM hydro
  subcycling, time-interpolated parent ghosts, per-child flux accumulation,
  reflux, and average-down at every branch;
- a four-level `1/2/3/2`-patch reactive conservation, synchronization,
  positivity, closure, and exact subcycle-count gate;
- symmetric chemistry--hydro--chemistry splitting over every patch-tree node,
  with deepest-to-root post-reaction synchronization and whole-tree rollback;
- a chemistry-on versus hydro-only branch comparison proving reactive species
  evolution while retaining composite mass, momentum, and total energy;
- recursive molecular transport with cumulative `r^2` child subcycling,
  time-interpolated parent ghosts, per-child diffusive flux registers, reflux,
  and average-down at every branch;
- full reaction--transport--hydro--transport--reaction composition with a
  four-level transport-on versus transport-off conservation gate;
- transactional runtime patch-tree rebuilds from supplied branching plans,
  with physical-coordinate overlap transfer across changed parent ownership;
- no-op, moved-tree conservation, deepest exact state/temperature retention,
  counter preservation, and invalid-plan rollback gates;
- per-parent normalized-gradient tagging and deterministic disconnected-tag
  clustering for automatic arbitrary-depth branching plans;
- tag-driven transactional tree rebuilds, including root-only creation,
  maximum-depth branching, unchanged-plan no-op, and invalid-request gates;
- independently owned adjacent children with parent-local same-level exchange
  for face and four-layer PPM/WENO ghost data;
- single owned time-integrated hydro and diffusive fluxes at fine/fine faces,
  with those internal faces excluded from coarse/fine reflux;
- adjacent PPM hydro and molecular-transport conservation, synchronization,
  exact exchange, and subcycle-accounting gates;
- deterministic cell-, hyperbolic-subcycle-, or parabolic-subcycle-weighted
  MPI ownership for every tree patch, collective hierarchy/work-model
  consensus rejection, owner-authoritative patch synchronization, and
  four-layer cross-rank adjacent-sibling halo gates;
- owner-only patch-tree chemistry with one global advance per patch, serial
  field parity, deepest-to-root synchronization, conservation, and global
  transactional rollback gates;
- owner-only recursive MPI patch-tree hydro with exact per-owner subcycle
  accounting, cross-owner adjacent PPM flux reconciliation, serial field
  parity, conservation, and global transactional rollback gates;
- owner-only recursive MPI patch-tree molecular transport with exact
  parabolic subcycle accounting, cross-owner shared diffusive fluxes, serial
  parity, conservation, and global transactional rollback gates;
- a transactional owner-only MPI `R-T-H-T-R` interval with complete
  bookkeeping synchronization, serial full-field parity, exact operator-call
  accounting, missing-database rejection, and rollback after a later-stage
  failure;
- rank-local sparse MPI AMR patch payloads, exact owner scatter/gather,
  one-copy global storage accounting, and same-hierarchy owner migration;
- direct owner-only chemistry on sparse AMR payloads with distributed
  average-down, parent/child ghost fill, adjacent PPM ghost replacement, and
  exact collective rollback;
- targeted point-to-point sparse MPI EB restriction from each child owner only
  to intersecting root tile owners, with exact transfer accounting and
  transactional rollback;
- owner-tiled direct sparse MPI EB hydro with point-to-point six-row halo
  exchange, bounded tile-local work, direct tile-to-child state/flux and
  correction routing, zero root-result/scatter traffic, exact accounting, and
  serial parity;
- owner-tiled direct sparse MPI EB SSPRK2 transport with point-to-point
  six-row halos, seam-isolated finite periodic-edge bands, targeted
  result/scatter and child traffic, zero-traffic tile-local final blending,
  exact accounting, and serial parity;
- patch-local EB flux-register storage and compact sparse MPI transport reflux
  correction round trips with unchanged transactional ordering;
- owner-local sparse MPI EB hydro/transport timestep selection with no root
  field traffic, fine-to-coarse subcycle scaling, communicator-minimum
  reduction, serial timestep parity, and collective rejection;
- public sparse MPI EB multi-step `R-T-H-T-R` advancement with a freshly
  selected stable interval per step, exact final-time clipping, committed-only
  clock and diagnostic publication, and serial full-field parity;
- transactional explicit sparse MPI EB topology rebuilding with serial
  overlap/prolongation parity, direct restriction/PCM/overlap owner traffic,
  deterministic owner recomputation, one-copy post-regrid storage, and
  complete invalid-control rollback;
- scheduled temperature-tagged sparse MPI EB topology rebuilding with
  root-owner planning, compact metadata broadcast, caller-defined geometry,
  direct owner migration, exact transfer accounting, serial dynamic-loop
  parity, and whole-step rollback;
- targeted root-only sparse MPI EB materialization with exact field parity,
  one packed send per remote root tile or child, unallocated non-root outputs,
  and collective invalid-payload rollback;
- writer-root-only sparse MPI EB formatted checkpoint and root/child CSV
  publication, serial checkpoint round-trip compatibility, collective I/O
  status, and exact successful-transfer accounting;
- root-only formatted checkpoint read with direct root-to-owner sparse restart
  scatter from a geometry-only replicated topology descriptor, exact
  per-remote-entity traffic, rank-local field parity, and collective
  metadata/I/O failure rollback;
- direct recursive hydro on sparse AMR payloads with mixed-ratio subcycling,
  replicated flux-register metadata, owner-local reflux/average-down,
  cross-owner PPM face reconciliation, and exact rollback;
- direct recursive molecular transport on sparse AMR payloads with cumulative
  `r²` subcycling, diffusive flux registers, cross-owner shared-face
  reconciliation, and exact rollback;
- direct SSPRK2 molecular transport on sparse MPI EB AMR payloads with
  owner-local child subcycling and reflux, distributed cut-interface
  conservation closure, serial parity, limiter parity, and exact rollback;
- a direct sparse `R-T-H-T-R` transaction with owner-only stage execution,
  exact call accounting, serial parity, missing-database rejection, and outer
  rollback after a later-stage failure;
- an end-to-end sparse MPI EB AMR `R-T-H-T-R` transaction that composes direct
  sparse chemistry, SSPRK2 transport, and hydro without a replicated child
  compatibility window;
- transactional topology-changing sparse regrid with rebuilt owner maps,
  exact overlap retention, one-copy persistent storage, serial parity, and
  invalid-plan rollback;
- tag-driven sparse regrid through four levels with disconnected-feature
  clustering, unchanged-plan no-op behavior, conservation, and invalid-tag
  rollback;
- packed point-to-point same-hierarchy owner migration with one direct message
  per changed patch and exact state, temperature, and ghost reconstruction;
- packed point-to-point adjacent sparse halo exchange with one bidirectional
  payload per cross-owner sibling face and no traffic on unrelated ranks;
- direct child-owner to parent-owner sparse interior transfer for chemistry
  average-down and hydro/transport synchronization;
- direct parent-state fanout only to distinct remote owners that require the
  state for sparse child ghost refresh;
- broadcast-free sparse recursive hydro and transport with owner-local flux
  registers, direct interval/flux/correction payloads, and reduced counters;
- an input-driven sparse MPI AMR driver with configurable subcycle-weighted
  ownership, stop-time clipping, periodic owner-local regridding, final
  conservation diagnostics, and ordered composite patch-tree CSV output;
- a moving-contact gate demonstrating lower AMR error than PCM.

For PCM/PLM, the reactive AMR application retains its overlap-preserving
two-level path when `amr_max_levels = 2`. PPM selects the multilevel engine at
any configured depth so every level owns the wider stencil state. The same
engine provides tag-driven arbitrary-depth state ownership, recursive
advancement, periodic hierarchy rebuilds, and composite output. A changed
multilevel hierarchy is conservatively averaged to the root before nested
patches are rebuilt, then old fine state and temperature data are copied
wherever old/new level spacing and physical cells align. Changed refinement
ratios fall back to conservative prolongation. With
`amr_multipatch_enabled = .true.`, the public application uses a tag-driven
two-level patch set, clusters disconnected tags, periodically rebuilds the
set, retains aligned fine overlap, and writes every uncovered parent or fine
cell exactly once. A separate static patch-tree engine permits each parent
patch to own zero or more separated children at arbitrary depth and advances
those patches recursively with level-ratio hydro subcycling. Every local child
accumulates its own coarse/fine boundary fluxes; reflux and average-down occur
after its fine subcycles, and failures restore the complete tree. The
patch-tree engine can synchronize to the root, tag and cluster every
prospective parent independently, and rebuild the resulting arbitrary-depth
branching plan transactionally. Independently owned adjacent siblings exchange
same-level ghosts, reconcile each shared time-integrated interface flux, and
exclude that internal side from reflux. The first MPI distribution bridge
assigns a unique owner to every patch and communicates authoritative fields
and adjacent halos while retaining replicas on all ranks. Chemistry and the
recursive hydro patch kernel execute only on those owners with collective
acceptance and rollback. Owner-authoritative start-state and face-flux
broadcasts let every replica apply the existing subcycling, shared-flux,
reflux, and average-down rules deterministically. Owner-only molecular
transport uses the same structure with `r^2` subcycling. One outer MPI
transaction now composes chemistry, transport, hydro, transport, and chemistry
and restores its synchronized starting tree after any rejected stage. The
sparse MPI bridge now removes non-owned patch payloads, reconstructs an exact
replica on demand, and moves payload ownership when a same-hierarchy owner map
changes. Chemistry is the first operator to run directly on sparse payloads:
child interiors are synchronized
deepest-to-root for average-down and parent/sibling data are streamed for
ghost refresh without materializing a complete tree. Recursive hydro now uses
the same sparse storage boundary while streaming interval states and fluxes;
molecular transport now follows the same structure with its parabolic `r²`
schedule. Their combined `R-T-H-T-R` transaction now remains sparse from its
outer backup through final acceptance, so normal full-physics advancement no
longer materializes a complete replica. Explicit-plan topology changes now
construct new patches on their owners, stream parent-prolongated interiors,
and transfer fine overlap directly from old owners to new owners. The same
sparse API now derives arbitrary-depth plans on parent owners, replicating only
compact integer topology metadata while building candidate fields on their
owners. Same-hierarchy owner changes now send one
packed patch directly from the old owner to the new owner. Adjacent sparse
siblings likewise exchange only the one- or four-layer boundary payload needed
by their two owners. Child interiors used by average-down and synchronization
now move directly from each child owner to its parent owner. Parent interval
start/end states now reach only distinct child owners, child boundary fluxes
return directly to the parent owner, and shared-face corrections return only to
the affected child owner. Level counters synchronize from owner-local deltas
once per physics stage. Final parent-to-child ghost refresh likewise sends one
parent state to each distinct remote child owner. Sparse physics and both
explicit-plan and tag-driven regrid contain no all-rank field replica. The
public sparse MPI AMR driver composes those APIs into a complete run; only the
initial root state and final diagnostic/output tree are materialized. A
periodic child may
touch a physical boundary only when it covers the full parent domain;
one-sided periodic refinement remains excluded because it crosses the periodic
seam.

## Build and test

Requirements: CMake 3.23 or newer, a Fortran 2018 compiler, and Python 3.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

MPI verification additionally requires an MPI implementation with the Fortran
2018 `mpi_f08` module:

```bash
cmake -S . -B build-mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_MPI=ON
cmake --build build-mpi --parallel
mpiexec -n 4 ./build-mpi/pelef_mpi_reactive_1d mpi_reactive_np4.csv
mpiexec -n 4 ./build-mpi/pelef_mpi_amr_reactive_1d \
  cases/mpi_sparse_amr_hotspot/hotspot.nml sparse_amr_np4.csv
mpiexec -n 4 ./build-mpi/pelef_mpi_eb_amr_patch_2d
```

To enable the live Cantera reference gate:

```bash
python3 -m pip install cantera==3.2.0
cmake -S . -B build-cantera \
  -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_CANTERA_REFERENCE=ON
cmake --build build-cantera --parallel
ctest --test-dir build-cantera --output-on-failure
```

With Ninja installed, the provided presets can be used for the ordinary suite:

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

## Representative runs

One-dimensional characteristic PLM:

```bash
./build/pelef cases/sod/sod_pelec_plm.nml
python3 tools/compare_sod.py --input sod_pelec_plm.csv
```

Two-dimensional periodic isentropic vortex:

```bash
./build/pelef2d cases/isentropic_vortex/vortex.nml
python3 tools/check_isentropic_vortex.py --input isentropic_vortex.csv
```

Passive two-species Sod problem:

```bash
./build/pelef_ms cases/multispec_sod/multispec_sod.nml
python3 tools/check_multispec_sod.py --input multispec_sod.csv
```

Elementary H2/O2 constant-volume reactor:

```bash
./build/pelef0d_h2o2 cases/zero_d_h2o2/reactor.nml
python3 tools/check_zero_d_h2o2.py --input zero_d_h2o2.csv
```

Reactive one-dimensional hotspot with characteristic PLM:

```bash
./build/pelef_reactive_1d cases/reactive_hotspot/hotspot.nml
python3 tools/check_reactive_hotspot.py --input reactive_hotspot.csv
```

The same case with monotone PPM and HLLC:

```bash
./build/pelef_reactive_1d cases/reactive_hotspot/hotspot_ppm.nml
python3 tools/check_reactive_hotspot.py --input reactive_hotspot_ppm.csv
```

The time-traced characteristic PPM path with the optional contact and shock
detectors enabled:

```bash
./build/pelef_reactive_1d \
  cases/reactive_hotspot/hotspot_characteristic_ppm.nml
python3 tools/check_reactive_hotspot.py \
  --input reactive_hotspot_characteristic_ppm.csv
```

Smooth general-EOS entropy wave:

```bash
./build/pelef_reactive_1d \
  cases/reactive_entropy_wave/entropy_wave.nml
```

General-EOS H2/N2 composition wave with HLLC:

```bash
./build/pelef_reactive_1d \
  cases/reactive_composition_wave/composition_wave.nml
python3 tools/check_reactive_composition_wave.py \
  --input reactive_composition_wave.csv
```


Dilute-gas transport coefficient probe:

```bash
./build/pelef_transport_probe transport_probe.csv
python3 tools/compare_transport_cantera.py --input transport_probe.csv
```

Periodic one-dimensional molecular-transport pulse:

```bash
./build/pelef_reactive_1d \
  cases/reactive_transport_1d/transport_pulse.nml
python3 tools/check_reactive_transport_1d.py \
  --input reactive_transport_pulse.csv --nx 96
```

Reactive two-dimensional hotspot with CTU and HLLC:

```bash
./build/pelef_reactive_2d cases/reactive_hotspot_2d/hotspot.nml
python3 tools/check_reactive_hotspot_2d.py \
  --input reactive_hotspot_2d.csv --nx 24 --ny 24
```

The same two-dimensional hotspot with characteristic PPM and the optional
contact/shock controls enabled:

```bash
./build/pelef_reactive_2d \
  cases/reactive_hotspot_2d/hotspot_characteristic_ppm.nml
python3 tools/check_reactive_hotspot_2d.py \
  --input reactive_hotspot_characteristic_ppm_2d.csv --nx 24 --ny 24
```

Oblique exact entropy-wave transport through the same 2D path:

```bash
./build/pelef_reactive_2d \
  cases/reactive_diagonal_wave_2d/diagonal_wave.nml
```

Oblique constant-pressure H2/N2 composition transport through the
characteristic-PPM/CTU path:

```bash
./build/pelef_reactive_2d \
  cases/reactive_diagonal_wave_2d/diagonal_composition_ppm.nml
```

Input-driven characteristic-PLM reactive EB hydro around a circular obstacle:

```bash
./build/pelef_reactive_eb_2d \
  cases/reactive_eb_circle_2d/uniform.nml
python3 tools/check_reactive_eb_circle_2d.py \
  --input reactive_eb_circle_2d.csv
```

Active-cell chemistry parity against the regular 2D path:

```bash
./build/pelef_reactive_eb_2d \
  cases/reactive_eb_chemistry_2d/reactive.nml
```

Runnable two-level reactive EB AMR hydrodynamics:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_2d/uniform.nml
python3 tools/check_reactive_eb_amr_2d.py \
  --coarse reactive_eb_amr_coarse_2d.csv \
  --fine reactive_eb_amr_fine_2d.csv
```

Two-level reactive EB AMR thermal conduction and diffusive reflux:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_transport_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_transport_2d/transport.nml
python3 tools/check_reactive_eb_amr_transport_2d.py \
  --reference-coarse reactive_eb_amr_transport_reference_coarse.csv \
  --reference-fine reactive_eb_amr_transport_reference_fine.csv \
  --transport-coarse reactive_eb_amr_transport_coarse.csv \
  --transport-fine reactive_eb_amr_transport_fine.csv
```

Temperature-tagged conservative fine-patch movement:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_2d/dynamic_hotspot.nml
python3 tools/check_reactive_eb_amr_dynamic_2d.py \
  --coarse reactive_eb_amr_dynamic_coarse_2d.csv \
  --fine reactive_eb_amr_dynamic_fine_2d.csv
```

Active-cell chemistry parity on both EB AMR levels:

```bash
./build/pelef_reactive_2d \
  cases/reactive_eb_amr_chemistry_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_chemistry_2d/amr.nml
python3 tools/check_reactive_eb_amr_chemistry_2d.py \
  --reference reactive_eb_amr_chemistry_reference_2d.csv \
  --coarse reactive_eb_amr_chemistry_coarse_2d.csv \
  --fine reactive_eb_amr_chemistry_fine_2d.csv
```

Static three-level reactive EB AMR:

```bash
./build/pelef_reactive_2d \
  cases/reactive_eb_amr_chemistry_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_2d/amr.nml
python3 tools/check_reactive_eb_amr_three_level_2d.py \
  --reference reactive_eb_amr_chemistry_reference_2d.csv \
  --root reactive_eb_amr_three_level_root_2d.csv \
  --middle reactive_eb_amr_three_level_middle_2d.csv \
  --finest reactive_eb_amr_three_level_finest_2d.csv
```

Static three-level checkpoint/restart parity:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_restart_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_restart_2d/checkpoint_stop.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_restart_2d/restart.nml
python3 tools/check_reactive_eb_amr_three_level_restart_2d.py \
  --checkpoint three_level_restart.chk \
  --reference three_level_restart_reference_root.csv \
    three_level_restart_reference_middle.csv \
    three_level_restart_reference_finest.csv \
  --stopped three_level_restart_stopped_root.csv \
    three_level_restart_stopped_middle.csv \
    three_level_restart_stopped_finest.csv \
  --restarted three_level_restart_restarted_root.csv \
    three_level_restart_restarted_middle.csv \
    three_level_restart_restarted_finest.csv
```

Tag-driven dynamic finest patch inside a fixed middle level:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_dynamic_2d/hotspot.nml
python3 tools/check_reactive_eb_amr_three_level_dynamic_2d.py \
  --root three_level_dynamic_root.csv \
  --middle three_level_dynamic_middle.csv \
  --finest three_level_dynamic_finest.csv
```

Reacting fine-to-root checkpoint/restart parity:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_restart_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_restart_2d/checkpoint_stop.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_restart_2d/restart.nml
python3 tools/check_reactive_eb_amr_restart_2d.py \
  --checkpoint reactive_eb_amr_restart.chk \
  --reference reactive_eb_amr_restart_reference_coarse_2d.csv \
  --stopped reactive_eb_amr_restart_stopped_coarse_2d.csv \
  --restarted reactive_eb_amr_restarted_coarse_2d.csv \
  --fine reactive_eb_amr_restart_reference_fine_2d.csv \
    reactive_eb_amr_restart_stopped_fine_2d.csv \
    reactive_eb_amr_restarted_fine_2d.csv
```

Reacting two-child checkpoint/restart parity:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_multipatch_restart_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_multipatch_restart_2d/checkpoint_stop.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_multipatch_restart_2d/restart.nml
ctest --test-dir build --output-on-failure \
  -R '^regression_reactive_eb_amr_multipatch_restart_2d_'
```

## Project records

- [Porting plan](docs/porting_plan.md)
- [Architecture](docs/architecture.md)
- [Numerical methods](docs/numerical_methods.md)
- [PeleC responsibility mapping](docs/pelec_mapping.md)
- [State variables](docs/state_variables.md)
- [Parity strategy](docs/parity_strategy.md)
- [Implementation status](docs/implementation_status.md)
- [Design decisions](docs/design_decisions/)


### PeleF 0.17.0 two-dimensional transport example

```bash
./build/pelef_reactive_2d cases/reactive_transport_2d/transport_hotspot.nml
python3 tools/check_reactive_hotspot_2d.py \
  --input reactive_transport_hotspot_2d.csv --nx 20 --ny 20
```


### PeleF 0.18.0 physical-boundary examples

```bash
./build/pelef_reactive_2d cases/reactive_boundaries_2d/couette.nml
./build/pelef_reactive_2d cases/reactive_boundaries_2d/thermal_channel.nml
./build/pelef_reactive_2d cases/reactive_boundaries_2d/inflow_outflow.nml
./build/pelef_reactive_2d \
  cases/reactive_boundaries_2d/prescribed_species_wall.nml
```

Solid walls are species-impermeable unless `wall_species_* = "prescribed"`
supplies a zero-sum species mass-flux vector in kg/(m2 s), positive from the
wall into the gas. Slip walls remove tangential viscous stress; no-slip walls
reflect velocity about a prescribed wall velocity. Prescribed fluxes provide
a catalytic-wall transport interface, not a surface-reaction-rate model.


### Full pressure-dependent H2/O2 chemistry

Set `chemistry_model = "full_h2o2"` to use the 10-species, 29-reaction third-body/Troe mechanism with the implicit cell reactor.


### PeleF 0.24.0 MPI 1D verification

Configure with `-DPELEF_ENABLE_MPI=ON`, then run the MPI verification drivers
with 1, 2, 4, or 8 ranks. The foundation and multispecies drivers use 257 cells
so the block decomposition is intentionally uneven; the chemistry, transport,
and coupled-reactive drivers use smaller non-divisible workloads to exercise the
same decomposition and ordered-gather logic.
