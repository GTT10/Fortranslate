# Regular 3D AMR patch tree

The existing `amr_patch_3d` describes a single child box in parent-cell indices.
`amr_patch_tree_3d_mod` composes those boxes into a tree instead of replacing the
qualified two-level kernels. Geometry/state operations are separate from time
integration so that multi-level conservation can be checked before adding flux
and subcycling logic.

## Representation and invariants

Node 1 is the root (parent 0, level 0). Every other node names an earlier parent;
`patches(i)` and `parents(i)` describe node `i+1`. Parent-before-child order is
required, but breadth-first order is not. Parent dimensions must match the
box's coarse dimensions. Siblings cannot share parent cells; face, edge and
corner adjacency are allowed. Refinement ratios are common at each depth and
may differ between depths. Physical domain boundaries are valid geometry.
There is no compile-time level or child-count limit.

The constructor checks bounds, parent relationships, overlap, level ratios and
integer extent overflow before replacing the destination. Each field stores a
four-dimensional `(component,x,y,z)` array. Components are generic conserved
quantities; this module does not infer thermodynamic temperatures or chemistry.

## Operations

- `initialize_amr_patch_tree_3d`: construct a topology transactionally.
- `prolong_amr_patch_tree_3d`: initialize all descendants with existing PCM
  prolongation from their parents; commit only when every node succeeds.
- `average_down_amr_patch_tree_3d`: traverse nodes in reverse order, reusing
  existing restriction so that descendants are synchronized before parents.
- `composite_integrals_amr_patch_tree_3d`: exclude cells covered by immediate
  children at every node, use the physical volume along its refinement chain,
  and accumulate all remaining cells with compensated summation. Covered
  parent values need not already be synchronized for this inventory.

Invalid topology/state leaves destination objects or inventory vectors unchanged.
The hierarchy owns its arrays; restriction uses a private candidate copy.

## Ghost cells and coarse-time interpolation

`fill_amr_patch_tree_ghosts_3d` builds a two-cell halo around one node. It
constructs checked 64-bit global integer coordinates at every depth, searches
all peer boxes at the requested resolution (including different parents), and
falls back through coarser levels when no peer covers the cell (PCM by default).
The result includes the unchanged interior and commits only on success.

Its component layout is density, three momenta, total energy, then any species
or other conserved components. Domain boundaries follow the existing regular
3D conventions: paired periodic faces, outflow, reservoir inflow, slip wall,
and adiabatic/isothermal no-slip walls. Euler mirroring flips normal momentum;
optional viscous mirroring includes tangential velocity and moving-wall work.
Corner face precedence matches the regular-grid helper. Temperature recovery
and isothermal heat flux are separate operations and are not performed here.

By default this operation samples one hierarchy snapshot. Optional paired
`coarse_future` and `coarse_fractions` arguments interpolate only donors below
the requested level. Fraction entry 1 belongs to the root, entry 2 to depth 1,
and so on; the vector length must equal the requested node depth. Every fraction
must be finite and in `[0,1]`. A donor uses its own depth's fraction, including
fallback past the immediate parent. Fractions 0 and 1 copy the exact endpoints.
All conserved components use the same convex weight. Boundary transforms are
applied after sampling/interpolation; reservoir inflow overrides donor values.

The future container must have one field slot per node, but only coarser slots
need allocated, shape-compatible, finite data. Peer/finer future fields are
unused and may remain unallocated. Current same-level peers and interiors are
always copied directly. Both optional arguments are required together, and
invalid interpolation input preserves the prior output.

Callers own time brackets and choose a fraction for each coarse level at the
fine stage time. They must supply synchronized coarse values where required.
Optional `coarse_limiter='mc'` or `'minmod'` enables conserved-component linear
reconstruction of coarse donors. Same-level copies remain exact. Slopes use
face neighbors, while extrema include the full 3x3x3 donor stencil to retain
multidimensional affine fields. Neighbor sampling crosses peer boxes and follows
physical/periodic boundaries, with PCM fallback if that stencil itself reaches
a coarser region. Temporal interpolation precedes slope construction and uses
each sampled donor's own time fraction.

One common slope scale bounds every component throughout the complete fine-cell
subdivision and preserves its zero-mean correction. Optional `species_count`
identifies rhoY components after energy; the final species slope is dependent on
density and the other species slopes, preserving their sum. Temperature, when
packed after the conserved fields, remains a guess for subsequent EOS recovery.
This component bound is not a proof of nonlinear EOS admissibility; existing
reactive flux/state checks can still reject a step. The recursive hydro API below
now supplies hydro subcycle scheduling; a separate transport scheduler uses ratio-squared substeps.
Coupled/application dispatch is now opt-in, as described below.
Peer lookup currently scans nodes; no scaling qualification is implied.

## Conservative flux synchronization

`reflux_amr_patch_tree_3d` accepts one `amr_tree_flux_3d` per node. Its x/y/z
arrays hold **time-integrated flux per unit face area over the same physical
interval**. Normal face indices are 1 through `n+1`; other indices are cell
indices 1 through `n`. The supplied state is the candidate already updated by
these unreconciled flux divergences. The routine does not multiply by a timestep.

Processing runs from deepest level to root. Peers first share one effective
face value, preferring a value already replaced by finer coverage. Equally
resolved peers must agree within `64*epsilon*max(1,maxabs(left),maxabs(right))`;
the earlier traversal value becomes canonical. This includes different-parent
peers and periodic seams. The divergence of effective minus original flux then
corrects uncovered cells. Child faces are area-averaged onto parents, including
faces on patch boundaries, so deeper corrections propagate upward. Covered
parent cells are supplied by the final bottom-up state restriction.

The returned `boundary_correction` is the signed physical-domain integral of
effective minus original root flux, to add to the original root boundary
ledger. Periodic faces contribute zero. `maximum_correction` measures the
largest uncovered-cell conserved correction, before average-down. Inputs are
read-only except for the candidate state and output diagnostics; all are
committed together only on success. Shape/nonfinite errors or conflicting
shared refined faces preserve prior state and diagnostics.

Global integer layout is shared with ghost filling. The implementation uses
node scans and candidate copies; scaling has not been measured. Flux production,
EOS admissibility after correction and subcycle scheduling remain caller work.

## Synchronous hydro

`advance_amr_patch_tree_hydro_3d` advances every node through SSPRK2 with a
common timestep. Each Euler stage fills snapshot halos, computes existing PCM
or characteristic-PLM Riemann fluxes, updates cells, reconciles integrated
fluxes, averages down and recovers temperature. The final convex combination
is restricted and recovered again. Initial covered parent cells are also
synchronized before flux evaluation. Final recovery uses the first-stage
temperature guesses, matching the regular-grid SSPRK2 convention.

State and one-component temperature trees, signed boundary change and maximum
reflux correction commit together only on success. The ledger averages both
Euler-stage surface integrals after fine-flux replacement. The maximum is
across the two Euler stages, before SSPRK averaging. Physical boundaries use
Euler reflection and optional reservoir state/temperature; thermal/viscous
wall transport is outside this routine.

`compute_amr_patch_tree_cfl_3d` returns the minimum existing local CFL bound over
all nodes, using each node's spacing. Callers choose the timestep; advance does
not enforce that bound. All nodes currently use this same interval. PLM face
reconstruction now also enables limited coarse-donor interpolation and forwards
the species count for closure. PCM reconstruction retains the original donor
path. Transport accepts the same optional coarse limiter. This is not a claim of globally second-order
AMR accuracy or subcycled evolution.

## Chemistry and reaction/hydro splitting

`advance_amr_patch_tree_chemistry_3d` applies the existing cell chemistry solver
to exposed cells at every node. Immediate-child boxes mask covered parent
cells; their result comes from bottom-up restriction. Temperatures are recovered
from the synchronized conserved state. Species and ordered reactions are
caller-supplied, with the existing optional chemistry-integrator selector.
The routine owns no fixed-mechanism database. State and temperature commit only
after all nodes and final recovery succeed.

`advance_amr_patch_tree_strang_3d` composes a half reaction interval, synchronous
hydro and a second half reaction interval. Candidate copies cover the whole
composition, so a hydro failure after chemistry also preserves the original
hierarchy and diagnostics. The returned boundary ledger is the hydro surface
integral; chemistry changes species internally, so conservation comparisons
use total mass, momentum, total energy and elements, not individual species.
The source-only interval may be zero; it still synchronizes covered parents and
recovers temperature. The coupled step requires positive time.

## Molecular transport

`advance_amr_patch_tree_transport_3d` performs common-timestep SSPRK2 transport.
Each stage packs conserved state and temperature guesses, fills viscous tree
halos and uses the existing physical ghosted flux kernel for viscosity, heat
conduction, species diffusion and optional barodiffusion. Domain-face wall
corrections are selected from global node coordinates; interior patch faces
receive no wall source. Reservoir, isothermal temperature and tangential wall
velocity inputs follow the existing regular-grid conventions.

Time-integrated fluxes pass through tree reflux before restriction and EOS
recovery. Final state, temperature, boundary ledger, largest reflux correction
and minimum species-flux limiter commit together. Shared-face incompatibility
or an inadmissible refluxed state rejects the step; automatic timestep retry
and an additional coarse/fine positivity limiter are not supplied here.

`compute_amr_patch_tree_transport_timestep_3d` returns the minimum local
transport bound and maximum diffusivity over all nodes. Optional domain faces
and wall temperatures also include isothermal-wall bounds on touched faces.
Callers must supply those options for isothermal walls and enforce the resulting
bound themselves. This remains a common timestep, not ratio-squared subcycling.

## Combined split update

`advance_amr_patch_tree_split_3d` composes R/2-T/2-H-T/2-R/2 using the tree
chemistry, transport and hydro APIs. Boolean switches skip chemistry or each
transport process; disabling all transport skips both transport stages.
Disabled operators do not require usable reaction/transport data. Barodiffusion
without species diffusion is rejected. Transport wall options are forwarded to
both half steps, while hydro uses Euler boundary conventions.

The boundary ledger sums both transport half-step ledgers and the hydro ledger.
Maximum reflux correction is the maximum of the three flow phases; minimum
species-flux theta is the minimum over the transport phases (one when disabled).
State, temperature and diagnostics commit only after the last reaction phase.
The caller still chooses a common timestep using the hydro, diffusion and wall
bounds and is responsible for retry after rejection.

## Conservative topology remapping

`remap_amr_patch_tree_3d` projects exposed old cells onto exposed new cells by
intersection volume. Root dimensions and physical extent remain the same;
patch locations, counts, depths and refinement ratios may change. Old covered
parent values are ignored, so they need not be synchronized. New covered parents
are filled by bottom-up restriction. This is piecewise-constant conservative
projection; it does not reconstruct subcell slopes.

Both topologies are represented on checked 64-bit least-common-multiple integer
coordinates. Face positions and intersection lengths are exact until the final
floating-point volume weights. Overflow is rejected before destination mutation.
Weighted components and coverage use compensated sums; each target cell must
have unit coverage within roundoff. The supplied old state/tree and new tree
are read-only, and the destination state commits only after restriction and
validation succeed. Distinct source and destination actual arguments are required.

`remap_amr_patch_tree_reactive_3d` additionally remaps temperature as an EOS
starting estimate, then recovers temperature from remapped conserved mixtures.
Both output containers commit together only on success. This API accepts the
new tree but does not replace the caller's active topology; callers publish
that topology only after success.

The implementation stores exposed old-cell descriptors and scans them for each
exposed new cell. Storage does not allocate a global finest-grid raster, but
work is quadratic in exposed cell counts. No scalable spatial search or timing
qualification is claimed.

## Serial application

`pelef_amr_tree_reactive_3d` accepts the existing flow namelist and an explicit
parent/box/ratio topology file. It samples initial conditions at each node's
global coordinates, then restricts and recovers parents. The driver chooses a
common all-node timestep including wall/inflow bounds and retries rejected
transactional steps by halving time. Exposed-cell CSV carries node/level/spacing;
a companion file records composite initial/final/boundary inventories. See the
[case guide](../../cases/amr_tree_reactive_3d/README.md).

The optional event schedule preloads valid topologies and clips timesteps to
their exact prescribed times. A candidate remap/EOS recovery and component
conservation check precede topology/state publication. The relative check uses
the old absolute-cell composite inventory, so opposite momenta do not create a
near-zero cancellation scale. Initial inventories and accumulated boundary flux
are preserved. Events at final time execute before output. Inflow bounds cover
the finest scheduled spacing; active-tree bounds are recomputed each step.

The configure-time selected frontend shares this driver through
`pelef_selected_amr_tree_reactive_3d_runtime`. Generated thermo/reactions/transport
and integrator choice are supplied by the frontend; the driver owns no fixed
mechanism lookup. Fixed callers omit the optional integrator and retain existing
automatic selection. Selected startup checks composition and the requested
initial/inflow/isothermal-wall temperature envelope against the generated data.

## Checkpoint storage API

`amr_tree_checkpoint_3d` groups topology, conserved/temperature fields, time,
step/regrid counters, initial/boundary inventories and reflux/limiter/remap
diagnostics. The versioned text writer validates finite state, parent-child
synchronization and EOS temperature consistency before opening the destination.
A physical write failure may leave a partial file; atomic replacement is not
claimed. The reader builds a private candidate and commits only after context,
geometry, data, footer and state validation all succeed. Stored doubles round-trip
without changing temperature guesses.

`build_amr_tree_checkpoint_context_3d` serializes flow configuration, ordered
NASA7/transport records and the existing complete reaction serializer, including
identity, composition, integrator and tolerances. Only final-time/step limits
and output path are excluded from configuration identity. The caller supplies
a canonical identity for initial topology and the full event schedule. Context
lines are compared exactly; formatting changes between compiler/runtime versions
may therefore reject a file even when its numeric inputs are equivalent.

The application now exposes optional `amr_tree_checkpoint` controls. Its context
includes canonical initial-tree geometry and each scheduled time/topology, without
filenames. Restarts validate run limits, the number of consumed events, expected
active topology and independently initialized inventories. The boundary ledger
uses `2e-10` times absolute cell inventories (a common vector magnitude for the
three momentum components); chemistry-enabled runs check the five Euler components,
while chemistry-disabled runs also check species. Successful-step checkpoints
are saved after same-time regrids and before final output. Exact saved fields,
temperatures, counters and diagnostics resume with the next unconsumed event.
The primitive API does not independently infer an event schedule from its opaque
caller-supplied identity or prove historical boundary-budget correctness.

## Temperature-driven topology generation

`tag_amr_patch_tree_temperature_3d` generates a new tree from exposed old-cell
temperatures. Callers supply nondecreasing positive temperature thresholds and
refinement ratios by depth, a parent-cell tile shape, a nonnegative tag buffer,
a node budget (including root), and periodic flags per axis. The number of
thresholds is the maximum refinement depth; an empty array requests root-only
geometry. A tile is selected if its buffered volume intersects an exposed cell
at or above that depth's threshold. Parent tiles partition the parent grid;
clipped final tiles handle shapes not divisible by the requested tile size.

Selection uses an exact integer common lattice for old and requested new levels.
It does not average temperature before tagging: a hot fine cell remains visible
even when the coarse parent is cold. Covered parent temperatures do not vote.
Buffered queries see cells across parent boundaries and wrap periodic axes;
physical domain faces do not wrap. Breadth-first parent order and k/j/i tile
order make the output independent of old patch ordering. Multiple separated
regions create separate children; absent tags remove old children.

The complete candidate is validated before publication. Invalid controls,
coordinate overflow, allocation failure and exhausted node budget leave the
caller's tree unchanged. This API creates geometry only; the existing conservative
reactive remapper transfers fields and recovers temperature. It scans hot cells
for each candidate tile, so large-tree scaling remains unqualified. Tile selection
is not a Berger-Rigoutsos clustering implementation or a proper-nesting guarantee.
The application now supports an optional `amr_tree_adaptation` group with
successful-step cadence. The first check follows a positive step count, preserving
initial-topology inventory semantics. Same-time updates precede checkpoint writes
and final output. Identical topology skips remapping but counts the check.
Prescribed schedules and automatic adaptation are mutually exclusive modes.

Checkpoint context includes the initial tree and all active adaptation controls.
On restart the counter must equal integer division of saved steps by cadence;
generated geometry must satisfy depth/ratio/tile/node bounds. A snapshot before
the first check must still match initial geometry. Retagging the saved temperature
is deliberately avoided: remapping can alter a formerly tagged cell's temperature.
The saved geometry and fields are restored and the next due step performs the
next check. This validates structural compatibility, not a replay of historical
tag decisions. Reservoir timestep bounds include the maximum requested depth.

## Transport coarse donors

`advance_amr_patch_tree_transport_3d` accepts optional `coarse_limiter`, forwarding
it and the species count to the tree ghost helper along with viscous-wall and
moving-wall data. The direct API defaults to PCM when the argument is absent.
The combined R/2-T/2-H-T/2-R/2 step supplies the configured limiter for both
transport half-steps when `reconstruction='characteristic_plm'`; PCM reconstruction
retains the previous path. Thus the application reconstruction selection now
controls spatial coarse donors consistently for hydro and molecular transport.

The existing transport flux routine recovers primitive/thermodynamic states from
the packed conserved ghost values and temperature guesses. Its positivity limiter,
physical wall fluxes, tree reflux and transactional commit are retained. Limiting
component interpolation alone does not guarantee nonlinear EOS admissibility;
failures still reject the step. This remains a common timestep across levels.

## Single-depth hydro stepping for subcycling

`advance_amr_patch_tree_hydro_level_3d` advances every patch at one requested
zero-based depth with synchronous SSPRK2 and returns time-integrated face fluxes.
The state, temperature and returned flux array commit together only on success.
Other depths remain unchanged, and their flux entries are unallocated. This
primitive deliberately performs no restriction or coarse/fine reflux; its caller
must accumulate fine substeps and synchronize their fluxes over a shared interval.

For a non-root depth, callers provide `coarse_start`/`coarse_end` endpoint trees,
packing conserved components followed by one temperature guess. Only coarser
fields must be allocated. `fractions_start` and `fractions_end` give the local
step endpoints within each donor level's interval; each array has length `depth`,
and each finite fraction pair must satisfy `0 <= start <= end <= 1`. Ghost values
use the start fraction for the first Euler flux and end fraction for the second.
Same-depth peers always come from the same current RK stage. PCM/PLM donor and
physical boundary behavior reuse the existing tree ghost/flux kernels.

The returned face arrays contain `dt/2 * (F_stage1 + F_stage2)` per unit area,
with the existing tree-register indexing. Their discrete divergence reproduces
the complete accepted RK2 update. The method recovers active temperatures after
the first stage and final convex combination; the first-stage temperature is the
final recovery guess. This API enables a recursive scheduler but does not itself
choose substep counts, enforce per-level CFL, reconcile levels or update the
composite boundary ledger. Until those pieces are connected, the application
continues to use its existing common-timestep path.

## Interval-scoped reflux

`reflux_amr_patch_tree_3d` now accepts optional zero-based `minimum_depth` and
an optional `synchronized_flux` output. The default depth zero retains full-tree
behavior. For a deeper minimum, only that depth and descendants are reconciled,
corrected and restricted; coarser states remain exact and their flux arrays may
be unallocated. The returned flux array contains effective, area-restricted peer
fluxes at active depths, with inactive entries unallocated. State, diagnostics
and optional flux output commit together only on success.

A nested subcycle can therefore reconcile its completed interval, accumulate
those effective fluxes over later substeps, then include its parent's original
flux when synchronizing the longer parent interval. Already corrected descendants
receive only any new effective-minus-supplied flux difference. Root boundary
correction is zero for a partial interval: its final physical boundary replacement
is accounted when the root interval closes. No asynchronous time intervals may
be mixed within one call; all supplied active fluxes must span the same interval.
This routine does not choose intervals, reconstruct temperatures or perform EOS
checks. Those remain the recursive hydro scheduler's responsibilities.

## Recursive hydro subcycling

`advance_amr_patch_tree_subcycle_hydro_3d` advances one coarse hydro interval
transactionally. Every depth performs synchronous peer RK2 steps. The next depth
performs `refinement_ratio` steps per parent interval, recursively; ratios 2 then
3 produce counts 1/2/6. Ancestor prediction endpoints remain available while finer
levels interpolate boundary values at both RK stage times. After a depth finishes,
interval-scoped reflux and restriction reconcile it with descendants, temperature
is recovered and effective fluxes return to its parent. Effective descendant fluxes
are summed over each parent interval. The final root boundary ledger combines its
original integrated surface flux and the fine-flux replacement.

Normalized interval endpoints share exact parent endpoints at the last substep.
Zero-width/underflowed intervals, failed RK stages, flux conflicts and failed EOS
recovery reject the entire coarse step. State, temperature, budget, maximum reflux
correction and optional `level_steps` counts commit together. The topology stays
fixed throughout the coarse interval. Scratch states and endpoints currently copy
full trees; this is not a memory/scaling optimization.

`compute_amr_patch_tree_subcycle_cfl_3d` returns the minimum local hydro CFL bound
multiplied by that depth's cumulative refinement ratio. The stepper accepts a
caller-chosen coarse `dt`; it does not retry or enforce this bound internally.
This hydro bound does not cover diffusion or external reservoir bounds.
The application and combined reaction/transport driver default to the common-step
path and can opt into the recursive composition. The transport scheduler below provides diffusive subcycling; coupled application
dispatch supplies the corresponding bounds and checkpoint controls.

## Recursive transport subcycling

`advance_amr_patch_tree_transport_level_3d` performs peer-synchronous transport
SSPRK2 at one depth, with viscous/thermal/species/barodiffusive fluxes, physical
wall data and coarse-time donor endpoints. It returns integrated RK fluxes and
the minimum positivity-limiter scale without restriction or reflux. Both Euler
states and the final average undergo EOS recovery. Final temperature recovery
uses the average of initial and second-Euler temperatures, matching the regular
transport convention. Nonactive depths and previous outputs remain unchanged.

`advance_amr_patch_tree_subcycle_transport_3d` combines that primitive with the
recursive endpoint/reflux sequence. A child executes `ratio**2` substeps, so ratios
2 then 3 produce counts 1/4/36. Integer counts and their products are checked for
overflow. The root physical budget includes fine boundary replacement; minimum
transport theta and maximum reflux correction accumulate over the full interval.
State, temperature, budgets and optional counts commit only on complete success.

`compute_amr_patch_tree_subcycle_transport_timestep_3d` scales local diffusive
and touched-isothermal-wall bounds by cumulative ratio squared, then takes their
minimum. Spatial spacing still follows cumulative ratio, not its square. The
reported maximum diffusivity is the maximum across nodes. As with hydro, the
caller supplies dt and remains responsible for stability, reservoir bounds and
retry policy. Full-tree scratch copies remain; this is not scaling qualification.
The combined driver and application can select these transport APIs through the
subcycling mode described below.

## Coupled application subcycling

The combined split accepts optional `subcycling` (default false). When true,
the two transport halves use recursive ratio-squared steps and hydro uses recursive
ratio-based steps. Chemistry retains the exposed-cell coarse half-interval solves.
All three flow budgets, maximum correction and minimum transport theta accumulate
inside the existing private candidate; the final chemistry solve must succeed
before publication. Disabled operators remain skipped.

The optional application `amr_tree_timestepping` group exposes this switch. Scaled
hydro/transport/wall timestep helpers replace the synchronous bounds, and reservoir
bounds use root spacing. The full dt is conservatively constrained by the transport
bound. Outer retries, scheduled-event clipping and automatic regrid cadence
continue to act on coarse intervals. Checkpoints occur only after full completion
and include the active integration mode in context; old synchronous contexts retain
their previous identity. Restart cannot change the mode. This is operator-split
subcycling, not per-fine-step coupled chemistry, and does not establish global
second-order convergence.

## MPI patch ownership and transfer

`mpi_amr_patch_tree_distribution_3d_mod` provides collective scatter/gather on a
caller-supplied communicator and root rank. Topology and an explicit zero-based
owner vector are replicated metadata. Each rank's local `amr_tree_state_3d` has
one field slot per node, but only owned nodes allocate cell arrays. Empty ranks
are supported. These partial fields intentionally do not satisfy the serial
full-tree `is_valid` contract and cannot be passed to serial evolution kernels.

Before payload messages, all ranks agree on transfer direction, root, node count,
exact topology/owner records and component counts. Shapes, finite data, ownership
allocation and MPI count bounds are checked collectively. Candidate allocations
also complete collectively before sends begin. Payloads travel directly between
root and each owner in node order. Scatter preserves its global source; gather
preserves local sources and changes the global output only on root. Failed input
validation leaves destinations unchanged on every rank.

The communicator is borrowed, and calls must follow collective ordering. Root and
owner buffers must be distinct actual objects. MPI runtime/communicator failure
is governed by the MPI error handler, not a recoverable application transaction.
These transfers support initialization/output and verification; they are not a
root-gather stepping algorithm. Ownership planning, evolution, rebalancing and application integration
remain open; distributed halos and restriction are described below. Whole-node
ownership does not split a large root patch or establish scaling.

## Sparse MPI cell requests

`fetch_mpi_amr_tree_cells_3d(tree,owners,local,nvar,requests,values,comm,ok)`
uses the same replicated metadata and owner-only payload layout. Each integer
request column is `[node,i,j,k]` with one-based indices. The component count is
common across ranks; request counts/order may differ, include duplicates, or be
zero. The committed output has shape `(nvar,number_of_requests)` in caller order.

After collective topology/ownership/shape/count checks, Alltoall exchanges request
counts and Alltoallv routes coordinates to owners. Owners return only those cell
vectors through a second Alltoallv. No full field gather is involved. Requested
values must be finite; unrequested cells are not scanned. All ranks validate and
allocate candidates before exchanges; invalid requests, metadata, allocation or
requested payloads preserve existing outputs. Communicator ordering and runtime
error handling follow the distribution API contract.

This primitive does not choose ghost donors or implement boundary transforms,
coarse spatial/time interpolation, restriction or reflux. The synchronous halo wrapper below constructs batched requests and consumes
the returned values; distributed evolution still requires further integration. Whole-node ownership and replicated metadata remain.

## Synchronous MPI halos

`fill_mpi_amr_tree_ghosts_3d` fills every owned node's padded `(nvar,nx+4,ny+4,nz+4)`
array, including its interior. Only owned output slots allocate. The API accepts
the same viscous-wall, inflow, wall-velocity, coarse-limiter and species-count
options as the serial fill, but currently has no coarse-time endpoint arguments.
Callers must supply consistent physical boundary options across ranks.

The serial implementation now delegates to `sample_amr_patch_tree_ghosts_3d`
with an `amr_tree_cell_source_3d` reader. This shares the actual geometry search,
periodic mapping, ancestor fallback, 27-cell PLM stencil, species-slope closure
and physical momentum/energy transforms. The serial wrapper retains full-state
and endpoint validation; its reader borrows source pointers only during the call.

The MPI reader first counts then records the donor sequence with planning enabled.
Planning walks geometry/stencils using finite placeholder values, skips slope and
wall arithmetic, and does not produce physical ghost values. A validation-only
mode checks options without allocating a root field on empty ranks. After the
sparse fetch, an evaluation walk consumes the cached sequence, verifying each
coordinate before use. All ranks agree on success before replacing output fields.
Both input rejection and post-exchange boundary-transform failure preserve outputs.

Requests currently include owned interior cells and repeated stencil donors;
there is no deduplication or persistent plan cache. Temporary storage scales with
local padded cells and requested stencil entries, rather than all tree payloads.
The planning walks and MPI self-requests are an optimization opportunity, not a
scaling result. Existing MPI runtime error-handler limits apply. The shared
sampler supports serial temporal endpoints, but this MPI wrapper is synchronous;
level-specific endpoint exchange and distributed evolution/reflux remain open.

## MPI bottom-up restriction

`average_down_mpi_amr_tree_3d(tree,owners,local,nvar,comm,ok,minimum_depth)`
updates covered parent cells without gathering full fields. The optional zero-based
minimum depth defaults to zero; only parents at that depth or deeper are changed.
The deepest level is a valid no-op limit. All ranks must agree on the depth.

The cell-exchange module now exposes its metadata/allocation preflight as
`validate_mpi_amr_tree_local_3d`; cell fetch and restriction both use it. Restriction
additionally checks all local values for finiteness, restricted-box MPI counts and
the shared kernel's default-integer `ratio**3` bound. This validation does not
claim that partial MPI states satisfy the serial full-tree validity contract.

Each rank privately copies only its owned payloads. Reverse node traversal updates
descendants before ancestors. The fine owner applies the existing
`restrict_average_3d` kernel to its current candidate and sends just the averaged
box to the parent owner (tag 28403). Same-owner pairs need no payload message;
other ranks allocate no box. Collective success gates each allocation, calculation
and transfer before the parent candidate is updated. Original local fields are
replaced only after the whole traversal succeeds.

This preserves the serial summation order and parent injection exactly, including
uncovered cells and parents above the selected depth. Calls borrow the communicator
and must obey collective ordering; MPI runtime failures follow its error handler.
The sequential node traversal uses per-node collective gates and is not a scaling
claim. EOS recovery, flux reflux, level stepping and application integration are
separate operations still to be connected.

## MPI peer flux synchronization

`synchronize_mpi_amr_tree_flux_peers_3d(tree,owners,nvar,depth,faces,flux,refined,comm,ok)`
reconciles shared faces at one zero-based depth. Active flux and `face_marks`
arrays are allocated only on each node's owner; other depths are untouched and
may be absent. The serial reflux module now exports its mark type. A true mark
means a value has already been replaced by finer-resolution flux.

The metadata preflight is factored separately from cell allocation checks so that
face operations need no dummy cell state. All ranks agree on topology, ownership,
component count, active depth and periodic axes; active local shapes, masks and
finite fluxes are validated before communication. The routine borrows its communicator.

Identical metadata derives intersections of peer boundary rectangles, including
periodic opposite faces and self-pairs. Each unordered node/face pair is processed
once. The higher-index owner sends values and encoded marks to the lower-index
owner (tag 28404); that owner reconciles them and returns the result (tag 28405).
Same-owner pairs bypass payload communication. If marks agree, the values must
match within the serial 64-epsilon tolerance and the lower-index/first-face value
is retained. Otherwise the fine-replaced value and mark propagate to both sides.
Physical boundary faces without peers and nonoverlapping portions remain unchanged.

Only owned active fluxes and marks are copied into candidates. Allocation,
validation or disagreement failures preserve all caller fluxes and marks; successful
completion moves active candidates back without touching other depths. MPI runtime
failures remain governed by the communicator error handler. Pair traversal and
collective gates are deliberately simple and have no scaling qualification.

This supplies the peer stage of distributed reflux. The composed API below connects
parent area averaging, uncovered-cell divergence correction and the boundary ledger.
The peer API alone does not update conserved state; EOS recovery remains separate.

## Distributed reflux composition

`restrict_mpi_amr_tree_flux_level_3d` transfers area-averaged child faces at one
positive child depth to parent owners. The child owner uses the serial transverse
summation order, dividing each contribution by the ratio twice. Each directional
averaged box travels directly to the parent (tag 28406). Parent candidates retain
existing fine-replaced marks; a second contribution must match the serial tolerance.
Only parent-depth fields commit, after every child/direction succeeds. Other depths,
including child source fluxes, are untouched.

`reflux_mpi_amr_tree_3d` composes peer synchronization and flux restriction deepest
first, applies effective-minus-original flux divergence to uncovered owned cells,
computes the root physical-boundary replacement ledger and calls distributed state
restriction. Its state input is the candidate already updated with the original
fluxes, just as for serial reflux. Input fluxes are time-integrated per face area.

The optional `minimum_depth` defaults to zero. Only that depth and deeper are
corrected/restricted; inactive input fluxes may be absent. Coarser state endpoints
remain unchanged and partial-interval boundary correction is zero. Optional
`synchronized_flux` returns owned active effective fluxes, with other slots absent.
The global boundary vector and maximum absolute cell correction are returned on
every rank using sum/max reductions. All ranks must agree on geometry, owners,
component count, floor depth, spacing and periodic geometry.

Only local owned states and active fluxes are copied. The wrapper keeps its own
candidates across all nested collective calls and commits state, effective fluxes
and diagnostics only when the entire correction succeeds. Thus a later-level peer
mismatch also rolls back earlier candidate work. Empty ranks participate in the
same collective sequence. No root gather is used for evolution or correction.

The public serial face access/match helpers and exported mark type avoid changing
serial tolerance or arithmetic semantics. Reflux checks finite states/fluxes but
does not recover temperature or prove EOS admissibility. Distributed RK stages,
coarse-time endpoint exchange and application dispatch must still be integrated.
Per-pair/per-child collective gates and whole-node ownership remain unqualified for
production scaling; MPI runtime failures follow the communicator error handler.

## MPI temperature recovery and common-step CFL

`recover_mpi_amr_tree_temperatures_3d` delegates every owned patch to the serial
reactive temperature-recovery kernel. An optional zero-based `depth` selects one
level; absent or -1 selects all levels. All ranks must agree on the selection.
Only selected owned temperature arrays are allocated as candidates, and successful
collective completion replaces them. Conserved input and unselected temperatures
are unchanged. Any owner recovery/allocation failure preserves all output guesses.

`compute_mpi_amr_tree_cfl_3d` computes each owned patch's serial CFL bound with
its recursively refined spacing and reduces the minimum across the communicator.
Empty ranks contribute the largest finite value. Spacing and CFL controls must
agree across ranks. Invalid inputs or EOS failure preserve the caller's dt.
This is the synchronous common-step bound, not a subcycle-scaled coarse bound.

Both APIs check replicated geometry/owners/component count, owner-only allocation,
finite conserved values and positive finite temperature guesses. Species records
must be identical on all ranks as a caller contract; component-count agreement
does not verify every thermodynamic coefficient. Partial-depth recovery still
validates structural and finite-value inputs at all depths. MPI runtime failures
follow the communicator error handler. These helpers do not advance time or
connect the MPI tree application.

## Synchronous MPI hydro stepping

`advance_mpi_amr_tree_hydro_3d` mirrors the serial common-step SSPRK2 algorithm on
owner-only state/temperature fields. The initial state is restricted and its
temperature recovered. Each Euler stage packs local conserved variables plus
temperature, fills distributed halos, evaluates the existing PCM or characteristic
PLM face kernels, updates local cells, applies distributed reflux, and recovers
temperature. The final half initial/half second-Euler combination is restricted
and recovered using the first-stage temperature guess, matching the serial order.

The root owner computes the original surface budget and broadcasts only that
component vector. Each stage adds the distributed reflux boundary replacement;
the returned RK budget averages both stages and the diagnostic takes their maximum
correction. Full fields are never gathered for stepping. Selected donors use the
same MC/minmod limiter and dependent species slope as serial PLM.

Topology/ownership/shapes and positive temperature guesses are checked collectively.
Ranks must agree on spacing, dt, boundary/reconstruction/limiter/solver controls and
active inflow arrays. Thermodynamic coefficient identity remains a caller contract.
Local allocations and flux-kernel failures are gated before the next communication.
All work stays private until both stages and final recovery succeed; failures preserve
state, temperature, budget and maximum correction. MPI runtime errors remain governed
by the communicator's error handler.

The caller selects dt (the common-step CFL helper is separate); there is no internal
retry or CFL enforcement. This API is hydro only, without chemistry, molecular
transport, level subcycling or application dispatch. Replicated metadata, duplicate
ghost requests and whole-node ownership retain the existing scaling limitations.

## MPI exposed-cell chemistry and Strang composition

`advance_mpi_amr_tree_chemistry_3d` invokes the serial reactive chemistry kernel
only on each owner's exposed cells. Direct-child boxes mask covered parent cells;
after every owner succeeds, distributed averaging and temperature recovery update
covered parents. Explicit/implicit selection is forwarded unchanged, including
the existing behavior when the optional integrator is absent.

Ranks agree on interval, tolerances, integrator selection and species/reaction
counts before calculation. Complete thermodynamic/reaction record identity remains
a caller contract. Owner-only shapes, finite conserved/temperature arrays and
positive temperature guesses are validated. Invalid input, allocation, ODE or EOS
failure is collective and leaves caller state/temperature unchanged.

`advance_mpi_amr_tree_strang_3d` holds private owner-local candidates across reaction
half-step, synchronous hydro step and reaction half-step. Hydro contributes the
boundary budget and maximum correction; chemistry preserves the existing total
energy/source convention. Outputs commit only after the final reaction succeeds,
so failure in hydro also discards the preceding successful chemistry half-step.
These are common-step operators without molecular transport, subcycling or MPI
tree application dispatch. MPI runtime failures retain the communicator's error
handler semantics; no full field gather or global chemistry raster is introduced.

## Synchronous MPI molecular transport

`advance_mpi_amr_tree_transport_3d` mirrors the serial transport SSPRK2 operator
on owner-local state and temperature fields. It supports viscosity, heat conduction,
species diffusion and pressure diffusion, with optional limited coarse donors.
Packed distributed ghosts use viscous reflection and moving-wall shifts. The
existing physical transport kernel sees only the physical faces touched by each
patch, derived from the integer tree layout. Other patch faces remain internal.

Each Euler stage updates owned cells, applies distributed reflux and recovers
temperature. The final convex combination is restricted; temperature recovery uses
half the initial plus half the second-Euler guess, matching serial transport.
Boundary/reflux budgets use the existing owner reductions, and minimum transport
limiter theta is reduced over every owner and both stages. Empty ranks contribute
one to this minimum.

Ranks must agree on geometry, dt/spacing, boundary faces, enabled terms, optional
wall/limiter presence, wall values and active inflow arrays. Transport and species
record identity remains the caller's responsibility; local kernel validation rejects
invalid records collectively. Only private candidates are changed until the entire
step succeeds, including final EOS recovery. State, temperature and all diagnostics
are preserved on rejection. There is no full-field gather.

This API advances transport alone with caller-supplied dt. A distributed transport
stability bound, full reaction/transport/hydro composition, temporal coarse donors
and application dispatch remain to be connected. Whole-node ownership, repeated
ghost planning and per-interface collective gates retain existing scaling limits.

## MPI transport bounds and full common-step composition

`compute_mpi_amr_tree_transport_timestep_3d` evaluates the serial diffusion bound
on each owned patch and the isothermal-wall bound only for the physical faces
that patch touches. It reduces dt by minimum and maximum diffusivity by maximum;
empty ranks contribute the largest finite dt and zero diffusivity. Geometry,
controls, enabled terms and optional boundary/wall inputs agree collectively.
On invalid input or local kernel failure both output scalars remain unchanged.
This is a common-step bound, without ratio-squared subcycle scaling.

`advance_mpi_amr_tree_split_3d` composes R/2-T/2-H-T/2-R/2. Disabled chemistry or
transport terms are skipped; enabled flags must agree before branching into
collective operators. Characteristic PLM selects the same limited coarse donor
for transport as for hydro. The wrapper holds local candidates across all stages,
sums the three flow/transport boundary budgets, takes maximum reflux correction
and minimum transport theta, and commits only after the final active stage succeeds.
A later hydro failure therefore discards preceding successful reaction/transport
work and leaves every caller diagnostic unchanged.

The chemistry module's existing validated owner-local clone helper is now shared
with the full composition; its numerical behavior is unchanged. Coefficient record
identity remains a caller contract. This API does not choose or enforce dt; the
caller combines hydro/transport/wall/inflow limits and handles retries. It exposes
synchronous operators, not recursive level subcycling, rebalancing or an MPI tree
application. Existing MPI runtime error-handler and scaling boundaries apply.

## Static fixed-mechanism MPI application

`pelef_mpi_amr_tree_reactive_3d <input.nml> <topology.dat>` uses a dedicated MPI
application driver. Metadata and configuration are replicated; whole node owners
are `mod(node-1,nranks)`. Analytic initialization samples each owned node's physical
coordinates directly, then distributed restriction and EOS recovery synchronize it.
Empty ranks still resolve reservoir composition and participate in collectives.

The driver combines global hydro/transport/wall bounds with reservoir bounds and
final-time clipping, then uses the common-step full split with up to 20 halvings
on failure. State is never gathered during evolution. Initial inventory and final
CSV/budget output deliberately gather to rank zero, preserving serial summation and
row order. This is a bounded output path, not distributed I/O or scalable partitioning.

Configuration/topology read errors and terminal step/I/O failures abort the supplied
communicator so another rank cannot remain waiting. The fixed frontend loads the
same elementary/full-H2O2 records as serial. The driver rejects tree extension
namelist headers rather than silently dropping checkpoint, adaptation or subcycling
requests. Schedule arguments are not accepted yet. Relative output/input alias
checks use POSIX getcwd before lexical normalization (1024-character capacity);
symlink/hard-link identity is not checked. Failed output writes can leave partial
CSVs, which are not restart files.

MPI tree checkpoint/restart, scheduled/automatic regridding, selected-runtime linkage
and recursive level stepping remain separate integration work.

## Remaining integration

MPI decomposition must be connected before claiming a general multi-level solver. EB geometry is separate.
The existing serial/MPI two-level application remains unchanged.

## MPI tree checkpoint continuation

The static fixed/selected MPI driver reuses the serial tree checkpoint schema and
canonical physics/mechanism/topology context. Rank zero gathers owned state and
stored temperature guesses only at requested checkpoints, writes through the
existing validator, and synchronizes completion. On restart rank zero validates
run limits, static geometry, initial inventory and boundary ledger in addition
to the reader's EOS/average checks, then scatters state and temperatures to the
current round-robin owners and broadcasts counters/diagnostics. No rank identity
is serialized. Temporary global payloads are released before stepping resumes.

The selected frontend passes its bundle SHA-256. Static serial/MPI interchange
therefore retains mechanism identity and exact temperature guesses. Scheduled or
automatic regridding and recursive modes remain outside this MPI entry point.
The existing writer's partial-file-on-I/O-error limitation also applies here.

## Sparse MPI tree remapping

`remap_mpi_amr_tree_3d` accepts old/new replicated tree metadata, independent
owner maps, an explicit field count and a borrowed communicator. It stores only
owned input/output payloads. Each rank builds the old exposed-cell geometry in
exact integer least-common-multiple coordinates, marks cells intersecting its
new exposed cells and requests each needed old cell once through the sparse
cell-query API. No full-state gather or replicated field payload is used.

Destination integration follows the serial source-cell order, overlap weights,
compensated accumulation and coverage normalization. Distributed bottom-up
restriction fills covered parent cells. A private candidate is published only
after all ranks pass planning, allocation, finite-value and restriction checks.
The reactive wrapper remaps conserved state and temperature guesses separately,
then recovers temperatures from the remapped mixture before committing either.
Source states remain unchanged. Species coefficient identity is a caller contract.

The serial `covered` and `gcd` geometry helpers are shared without changing serial
arithmetic. Integer-coordinate overflow is rejected. This initial implementation
replicates exposed-cell geometry and scans old/new cell intersections quadratically;
it does not claim production scaling or cache remap plans. Application event/tag
selection, geometry replacement, owner-map policy and checkpoint event metadata
are a subsequent integration layer.

## Scheduled MPI application regridding

The fixed/selected MPI application accepts the serial schedule file as an optional
third argument. It preloads all event geometries, includes the complete schedule
in checkpoint context, and bounds reservoir timesteps for the finest scheduled
spacing. At each event the sparse reactive remapper transfers old exposed cells
to round-robin new owners. Conservation uses serial-order compensated inventories:
each node owner continues the accumulated sum/correction and broadcasts only those
vectors. The same 5e-12 scaled remap threshold applies without full-field gathers.

Events run before checkpoint/output at their exact clipped times. Restart validates
the consumed event count and expected event geometry, then reconstructs metadata
on every rank before scattering checkpoint state to current owners. It restores
remap diagnostics and resumes at the next event without replay. Serial/MPI formats
and mechanism identities remain shared. Temperature-driven event generation and
recursive timestepping are still separate, unfinished MPI integration work.

## Distributed temperature tagging

`tag_mpi_amr_tree_temperature_3d` validates replicated geometry/owners and exact
agreement of thresholds, ratios, tile extents, buffer, node limit and periodic
axes. Each owner identifies hot exposed cells on the same checked integer lattice
as serial. Only hot-cell boxes and temperatures are all-gathered; cold payload and
covered parent cells do not vote. The extracted `build_tagged_tree_3d` helper
constructs the same breadth-first k/j/i tile order for serial and MPI. It is an
internal shared builder whose callers validate its controls/lattice.

Source order in the gathered hot-cell list is irrelevant because each tile takes
a Boolean intersection vote. Empty owners and zero requested refinement levels
are supported. Input/coordinate/allocation/node-limit failures are collective;
only an all-rank successful candidate replaces the destination tree. Temperature
source payloads are read-only. MPI counts are bounded before packing six integer
box coordinates per hot cell.

Hot-cell geometry/temperature is replicated, so this is not a production scaling
claim for nearly all-hot domains. Construction still scans tiles against hot
cells. Automatic application invocation and restart adaptation counters remain a
subsequent integration layer.

## Automatic MPI application adaptation

The fixed/selected MPI driver now accepts the serial adaptation namelist and
context records. Scheduled and automatic modes remain mutually exclusive.
Initial proposal validation checks lattice/node limits; the finest requested
spacing supplies conservative reservoir bounds. After each configured successful
step interval, distributed temperature tagging proposes geometry and the sparse
reactive remapper applies changes. The adaptation counter advances even when
geometry is unchanged, matching serial checkpoint semantics.

Adaptation precedes checkpoint/output, including the final step. Restart checks
`steps / regrid_interval_steps`, permitted depth/ratios/node count and tile bounds.
Only rank zero reads the file; validated parent/box/ratio metadata is broadcast
and rebuilt on every rank before owner-local state/temperature scattering. The
last completed step is recorded to prevent duplicate tagging at continuation.
Serial/MPI checkpoint context and payloads remain compatible for common timesteps.

## MPI coarse-time ghost donors

`fill_mpi_amr_tree_ghosts_3d` now accepts paired `coarse_future` and
`coarse_fractions` arguments, plus optional zero-based `depth`. Omitted depth or
`-1` fills all owned nodes; a nonnegative depth allocates output only for owned
nodes at that depth. Fractions contain one value per coarser depth, or
`tree%level_count()-1` entries for all-depth filling. Each sampled node receives
the prefix matching its depth. Presence, depth and finite [0,1] fractions must
agree across ranks. Same-level peers continue to use current state.

Planning records both cell coordinates and whether each read is current/future.
Separate sparse exchanges populate those requests, then the serial sampler applies
PCM/limited PLM time blending and physical-boundary transforms in the same order.
The current state is complete on owned nodes; future fields need only exist on
owned coarser levels. All required coarse future values must be finite even when
fractions are zero, matching serial validation. Inactive finer endpoint fields
are ignored. Outputs commit only after all ranks finish sampling successfully.

The sparse cell fetch/local validator has an optional `maximum_depth` bound.
It validates only fields at or below that depth and rejects requests beyond it;
`-1` permits an empty payload/query set. All ranks agree on the effective bound.
Default calls retain full-tree validation. This permits partial future endpoints
without placeholder allocation or full-field communication. Numerical level
stepping and recursive flux accumulation still need to consume these donors.

## MPI single-depth hydro step

`advance_mpi_amr_tree_hydro_level_3d` advances only owned nodes at a selected
zero-based depth and returns owner-local time-integrated face fluxes. It reuses
the existing PCM/characteristic-PLM flux kernels, sparse temporal ghost sampler
and owner-local EOS recovery. Endpoints pack conserved components plus a stored
temperature guess; they are required for depths above zero and only coarser
owned fields must be allocated. Stage one/two use their respective fraction
vectors and retain current same-level peers.

The API follows serial SSPRK2 order: Euler stage one, active temperature recovery,
Euler stage two, half-initial/half-stage-two combination, then active temperature
recovery from the first-stage guess. Face integrals use half-dt times each stage
flux in serial arithmetic order. No averaging down, peer synchronization or
reflux occurs here; those belong to the enclosing completed level interval.
Inactive state/temperature fields are preserved and inactive flux slots remain
unallocated. All result objects commit only after successful collective completion.

Depth, dt/spacing, options and endpoint presence are agreed across ranks. The
shared ghost sampler checks stage fractions. A disagreement discovered at stage
two rolls back the complete call. Species coefficient identity remains a caller
contract. The caller supplies a stable dt; the API does not compute/enforce CFL
or retry failed steps. Recursive level scheduling remains separate work.

## MPI single-depth transport step

`advance_mpi_amr_tree_transport_level_3d` advances one zero-based depth using the
shared molecular-transport kernel, sparse coarse-time ghosts and owner-local EOS
recovery. Current/endpoint packing and controls match the hydro level API. It
supports viscosity, conduction, species diffusion and barodiffusion, optional
PCM/MC/minmod coarse donors, reservoir inflow and prescribed moving/isothermal
walls on only the physical faces touched by an active patch.

The serial transport ordering is retained: recover temperature after both Euler
stages, combine the initial and second-stage conserved states, average initial
and second-stage temperature guesses, then recover temperature once more. Face
flux integrals sum half-dt contributions from both stages. The reported minimum
positivity factor is the minimum across active owners and both stages; empty
owners contribute one. Inactive states/temperatures remain unchanged and their
flux slots are unallocated.

Private candidates cover the whole call, including late second-stage failures;
state, temperature, flux and minimum factor publish together. The API performs
no reflux, averaging down, CFL enforcement or dt retry. Recursive ratio-squared
transport intervals and interval flux synchronization remain to be connected.
Species/transport coefficient identity across ranks remains a caller contract.

## MPI recursive hydro intervals

`advance_mpi_amr_tree_subcycle_hydro_3d` now follows the serial recursive interval
schedule. It snapshots old/new packed endpoints on the owner of each advanced
node, computes ancestor fractions from normalized interval times, performs one
MPI level step and then advances the next depth by the refinement ratio. Child
integrated fluxes accumulate on their node owners in substep order. Each completed
interval runs depth-limited distributed reflux and temperature recovery before
returning synchronized effective fluxes to its parent interval.

Only the root owner integrates physical root face budgets; that small vector is
broadcast before the final root reflux correction. State, temperature, boundary
budget, maximum correction and optional per-depth counters publish together after
the complete coarse interval succeeds. Temporary endpoints and face fluxes stay
owner-local. No global field gather is introduced.

`compute_mpi_amr_tree_subcycle_cfl_3d` computes owner-local CFL bounds using
spacing divided by the cumulative refinement ratio and scales local dt back to
the coarse interval before a global minimum. Empty owners contribute huge dt.
Counts use checked int64 products. The advance API does not enforce this bound;
the caller selects dt/retry policy. Runtime coefficient identity remains a caller
contract. Recursive transport and application dispatch remain separate work.


## MPI recursive transport intervals

`advance_mpi_amr_tree_subcycle_transport_3d` follows the serial ratio-squared
schedule and transport level RK2 arithmetic. It retains owner-local packed
coarse endpoints, accumulates child face integrals in substep order and runs
depth-limited distributed reflux and EOS recovery at each completed interval.
The global minimum transport limiter joins the boundary budget, maximum reflux
correction, state, temperature and counters in the whole-interval transaction.
No field gather is introduced. Refinement-count products use checked int64.

`compute_mpi_amr_tree_subcycle_transport_timestep_3d` evaluates owner-local
transport and touching-isothermal-wall bounds at serial integer-layout spacings,
rescales by cumulative squared ratios, then reduces the minimum timestep and
maximum diffusivity. Controls, optional wall inputs and face choices must agree
across ranks. Empty owners contribute huge dt and zero diffusivity. Coefficient
identity remains the caller's contract; advance does not enforce dt or retry.
Recursive split composition and application mode dispatch remain to be connected.


## MPI recursive split composition

`advance_mpi_amr_tree_split_3d` accepts optional `subcycling`, default false.
Its effective value joins the collectively agreed operator flags before any
stage runs, preventing ranks from entering different communication schedules.
When enabled, both transport half-intervals use ratio-squared recursive steps
and the full hydro interval uses ratio-based recursive steps. Chemistry retains
the existing exposed-cell coarse half-intervals. The order is R/2 T/2 H T/2 R/2,
matching the serial split API. Disabled operators are skipped, and PLM forwards
the configured coarse-donor limiter to transport.

The outer candidate transaction retains all state and diagnostics until every
stage succeeds. Failure after earlier chemistry or transport work discards that
work. Default/false mode retains common-step arithmetic. Application namelist,
scaled bound selection and checkpoint mode identity are the next integration.


## MPI application recursive mode

The fixed/selected MPI application reads `amr_tree_timestepping` with the same
`subcycling` default and semantics as serial. It selects scaled hydro and
transport/wall bounds and forwards the mode into the complete split transaction.
Reservoir bounds use root spacing in recursive mode, since substep counts cancel
the spatial refinement factors. Common-step bounds remain unchanged. Retry,
event clipping, adaptation checks and checkpoint cadence count complete coarse
intervals; no unfinished substep is published or saved.

The canonical schedule/adaptation identity appends `subcycling-v1` only when
enabled, at the same location as serial. Thus old common-mode contexts remain
compatible, recursive checkpoints can cross serial/MPI frontends and rank counts,
and a different mode is rejected before payload publication. Rank zero logs the
selected mode. Unknown tree groups and malformed timestepping inputs fail explicitly.
The focused application test covers a static two-level tree; dynamic recursive
application qualification remains separate.


## Work-weighted MPI whole-patch ownership

`partition_amr_tree_3d` assigns decreasing estimated patch work to the rank with
the smallest accumulated work. Equal weights prefer the lower node index and
equal rank loads prefer the lower rank. Work is cell count times cumulative
refinement ratio to exponent 0, 1 or 2 for common, recursive hydro or recursive
transport stepping. Counts, products and rank sums are checked int64 operations;
invalid input, allocation failure and overflow leave the prior owner vector intact.
No field payload is inspected. Replicated input geometry and controls are a caller
contract; existing transfer/step routines collectively validate ownership.

The MPI application accepts `amr_tree_distribution partitioning='weighted'`;
default `round_robin` preserves prior ordering. A single assignment helper serves
initialization, remap destinations and restored geometry. The policy is excluded
from physical checkpoint identity, permitting policy and rank changes on restart.
Existing sparse remap migrates data when geometry changes; no independent measured
load-triggered migration is introduced. The estimate excludes chemistry stiffness,
halo/communication costs and covered-cell differences. Whole-patch indivisibility
can dominate imbalance. Selection currently costs O(nodes squared + nodes*ranks),
with no large-tree scaling claim; spatial locality is not optimized.


## Owner-local final CSV

`write_mpi_amr_tree_csv_3d` uses the existing row formatter and exposed-cell
predicate. The first owner-local pass validates every output row and measures
per-node byte counts without opening the destination. A small metadata reduction
and checked MPI-offset prefix sum give each node its serial-order byte interval.
The second pass regenerates rows into bounded chunks and uses independent
`MPI_File_write_at` calls. The collective file size replaces any older tail;
rank zero writes the shared header. Empty owners participate in collective
validation/open/close without field payloads. Path, header, time and spatial
controls must agree across ranks. Short writes are rejected through status counts.

The application keeps final fields on their owners. Its existing compensated
inventory traversal passes only sum/correction vectors in serial node/cell order;
rank zero writes the small budget file. Initial inventory and checkpoint gathers
remain separate work. Preflight failures preserve the destination, while actual
I/O failure after open/truncation can leave a partial file; no rename transaction
or parallel-filesystem performance claim is made. Formatting twice trades CPU
work for bounded output memory and input validation before destructive I/O.


## Owner-local checkpoint saving

The MPI checkpoint writer accepts replicated tree/ledger/counter metadata plus
owner-local state and temperature. Before opening the destination it clones local
fields, performs distributed averaging and compares the result with input using
the serial 64-epsilon parent tolerance. The shared checkpoint EOS/temperature
consistency helper validates each owned patch. Metadata finiteness/counters,
context structure, identical path/header and checked byte offsets are also
validated collectively. Failure preserves the prior destination and input fields.

The header retains the V1 serial format; fixed-width state/temperature rows have
precomputed node intervals in the file. Owners write bounded chunks independently
with MPI-IO; rank zero writes metadata and the footer. The collective file-size
operation removes stale tails. No state gather or serial-writer fallback is used.
Temporary validation copies stay on their owners. Actual I/O failure can leave a
partial destination; this is not an atomic rename transaction. Loading still uses
the validated serial root reader and owner scatter. Checkpoint bytes and physical
identity do not depend on ranks or assignment policy.


## Indexed owner-local checkpoint loading

The serial metadata parser is shared through `read_amr_tree_checkpoint_header_3d`.
The MPI reader opens formatted streams and reads metadata on each rank. V1 stores
no byte index, so rank zero reads one list-directed cell record into a scratch
vector at a time, records each node's stream position and checks footer/EOF.
Only node positions are broadcast. This preserves the existing text/list-directed
record parsing rather than assuming every externally supplied record has the
writer's fixed byte width. It requires a serial scan and a second read of payloads,
but never retains a complete root state.

Each owner seeks to its node offsets and allocates/reads only its assigned U/T.
The same collective snapshot preflight as saving checks metadata/context, EOS
and parent averages. Owner arrays and snapshot outputs publish only after success.
The application then checks expected geometry, counters and initial/boundary
ledgers using distributed compensated inventories before replacing active fields.
The old root snapshot scatter and separate adaptive-geometry broadcast are removed.
Current partition policy determines restored owners; rank and policy changes remain
compatible with serial V1 checkpoints. Removing the index scan requires additional
indexed-format work and is not claimed by this memory-bounded reader.


## Optional indexed checkpoint V2

`checkpoint_indexed` selects V2 when saving; default false retains V1 bytes.
The V2 magic precedes the unchanged context/ledger/topology records, followed by
one int64 relative byte offset per node plus the payload-end offset. Offsets are
relative to the first cell record, avoiding dependence on header lengths. A
shared checked geometry/record-width helper builds the index for both writers
and validates it in the common header reader. Fixed-width records contain nvar+1
ES27.18E3 fields separated by one space and terminated with LF, so each cell
occupies exactly 28*(nvar+1) bytes.

Readers auto-detect the format. MPI adds the payload-start stream position to
the verified offsets and seeks to the footer directly; each owner then reads
its node ranges. Serial and MPI readers verify consumed bytes against node
intervals. Wrong offsets, overflow, wrong footer/EOF and noncanonical V2 record
lengths reject candidates. V1 keeps its existing list-directed compatibility and
streamed indexing path. Index presence is not a physical restart-context change;
new readers accept either format, but V1-only older executables cannot read V2.


## Direct same-geometry owner migration

`redistribute_mpi_amr_tree_3d` accepts immutable old/new owner maps and publishes
both conserved fields and temperatures together. Collective metadata/layout,
finite-state, positive-temperature and MPI count checks precede communication.
All destination payloads are allocated before the first send; allocation failure
is agreed collectively. Owners traverse nodes in the same order, copying unchanged
owners and sending changed owners directly with separate U/T tags. The candidate
replaces caller fields only after collective success. Nonowners retain empty
payload slots. No interpolation, EOS solve, restriction or root gather changes
values; the caller adopts the new owner map after success.

This primitive supports future runtime rebalancing but does not schedule it or
measure work. It temporarily holds old and candidate owner-local fields and uses
blocking, ordered whole-patch transfers. MPI runtime failures follow the caller's
error handler; process-loss recovery and production scalability are not claimed.


## Cost-driven MPI rebalance decision

`rebalance_mpi_amr_tree_3d` receives patch costs for a common completed sampling
interval, plus a replicated minimum fractional gain in [0,1). It only inspects
entries owned by the caller rank, requiring finite nonnegative values. An allreduce
collects these uniquely contributed costs. Dividing by the largest cost before
rank accumulation avoids overflow for large finite timings; all-zero cost is a
successful no-op. Largest-cost-first placement uses lower node/rank indices for
ties; zero-cost nodes retain their owners. The candidate is accepted only when
its maximum rank cost is strictly lower and the fractional decrease reaches the
threshold. Otherwise owners and fields remain unchanged.

The accepted plan calls direct U/T migration, then publishes the owner vector.
Invalid metadata, cost/threshold disagreement, allocation or migration preflight
failure leaves caller owners and fields intact. Shape validation occurs even on
no-op calls; thermodynamic consistency remains the caller's contract. This API
neither collects timers nor schedules runtime calls. It assumes comparable costs
across ranks and does not account for migration cost, topology, communication,
rank speed differences or sampling noise. Predicted load reduction is not measured
speedup. Planning retains quadratic patch selection complexity.
