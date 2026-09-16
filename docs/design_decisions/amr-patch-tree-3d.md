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

## Remaining integration

MPI decomposition must be connected before claiming a general multi-level solver. EB geometry is separate.
The existing serial/MPI two-level application remains unchanged.
