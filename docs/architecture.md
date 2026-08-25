# PeleF architecture

## Executable split

PeleF exposes twelve serial verification drivers and seven optional MPI drivers
over shared numerical and physical-property modules.

```text
pelef
  └─ one-dimensional constant-gamma PCM / PLM / PeleC-style paths

pelef2d
  └─ two-dimensional periodic constant-gamma CTU-style Euler path

pelef_ms
  └─ passive runtime-multispecies transport over the constant-gamma core

pelef0d
  └─ NASA7 mixture thermodynamics and a synthetic isomerization reactor

pelef0d_h2o2
  └─ generated reversible elementary kinetics and an H2/O2 reactor

pelef0d_h2o2_full
  └─ full pressure-dependent H2/O2 kinetics and an implicit reactor

pelef_transport_probe
  └─ qualified dilute-gas mixture transport coefficients

pelef_reactive_1d
  └─ NASA7 reactive Euler with PLM/PPM, HLLC/PeleC fluxes, and Strang splitting

pelef_reactive_2d
  └─ NASA7 reactive Euler, physical boundaries, transport, and CTU correction

pelef_reactive_eb_2d
  └─ reactive embedded-boundary hydro with weighted StateRedist

pelef_reactive_eb_amr_2d
  └─ lifecycle-capable reactive EB AMR hydro, chemistry, and restart

pelef_amr_reactive_1d
  └─ dynamic two-level reactive AMR with PLM, chemistry, and transport
```

The established constant-`gamma` solvers remain intact as regression baselines. Composition-dependent flow is introduced through a separate driver and modules so a general-EOS change cannot silently alter earlier results.

## Shared thermodynamics and chemistry

```text
thermo_database_mod
  └─ named NASA7 species records and molecular weights
        ↓
nasa7_thermo_mod
  └─ species cp, cv, h, u, and standard-state entropy
        ↓
mixture_thermo_mod
  ├─ mixture molecular weight and gas constant
  ├─ mixture cp, cv, gamma, h, and u
  ├─ ideal-gas pressure, density, and frozen sound speed
  └─ safeguarded e(Y,T) -> T inversion

mechanisms/h2o2_elementary.json
        ↓
tools/generate_elementary_mechanism.py
        ↓
src/generated/h2o2_elementary_mechanism_mod.F90
        ↓
elementary_kinetics_mod
        ↓
constant_volume_reactor_mod
```

The normalized JSON file is the authoring format for the current generated reaction subset. CI regenerates the committed Fortran module and requires byte-for-byte equality.


## Shared molecular transport

```text
transport_database_mod
  └─ pinned Lennard-Jones records for the seven-species subset
        ↓
mixture_transport_mod
  ├─ Chapman--Enskog pure viscosity
  ├─ Wilke mixture viscosity
  ├─ modified-Eucken / Mathur thermal conductivity
  ├─ Chapman--Enskog binary diffusion
  └─ mixture-averaged species diffusion
        ↓
reactive_diffusive_flux_x
  ├─ Newtonian viscous stress and viscous work
  ├─ Fourier heat flux
  ├─ mole-fraction and optional pressure diffusion driving forces
  ├─ correction velocity enforcing sum(j_k) = 0
  └─ species-enthalpy diffusion energy flux
```

Transport data are pinned to the same Cantera `h2o2.yaml` provenance as the
current thermodynamics/chemistry subset. The implemented coefficient model is
a dilute ideal-gas subset, not the full PelePhysics generated polynomial
transport layer.

## Reactive one-dimensional path

The reactive state is stored as

```text
state(variable, 0:nx+1)
```

with the conserved layout

```text
rho, rho*u, rho*v, rho*w, rho*E, rho*Y_1 ... rho*Y_N.
```

Temperature is a synchronized derived field, not an independently evolved conserved variable.

```text
state + temperature guess
        ↓
reactive_conserved_to_primitive
  ├─ recover Y from rho*Y
  ├─ remove kinetic energy from rho*E
  ├─ solve u(Y,T) = e_target
  ├─ evaluate p(Y,rho,T)
  └─ evaluate frozen sound speed
        ↓
reconstruction / Riemann / CFL
```

Hydrodynamic responsibilities are separated as follows:

```text
reactive_primitive_to_conserved
reactive_conserved_to_primitive
        ↓
reconstruction selector
  ├─ PCM
  ├─ characteristic PLM + MUSCL-Hancock tracing
  ├─ componentwise monotone PPM + SSPRK3
  └─ characteristic PPM profile integration
       ├─ optional contact steepening
       └─ optional shock flattening
        ↓
general-EOS Rusanov, HLLC, or PeleC-style acoustic flux
        ↓
conservative finite-volume update

optional transport branch
  ├─ face-centered viscous / conductive / species fluxes
  ├─ explicit SSPRK2 diffusion update
  └─ parabolic dx^2 / diffusivity timestep gate
```

Advective species face fluxes close to the total mass flux. Diffusive species
fluxes use a correction velocity so their sum is zero to roundoff.

The reactive `pelec` selection follows the upstream acoustic construction.
Each side supplies its NASA7 frozen sound speed and acoustic impedance. The
solver estimates star pressure and normal velocity, selects an upwind vector
of species densities (or averages it at a stationary interface), applies the
star-density pressure correction species by species, and evaluates the star
sound speed through the mixture EOS. Inward/outward wave interpolation selects
the final interface density, composition, velocity, and pressure. Total energy
is then rebuilt from that interface through the NASA7 EOS before assembling
the conservative flux. Rusanov and HLLC remain separate selectable kernels.

## Reactive two-dimensional CTU path

The reactive 2D state is stored as `state(variable,nx,ny)` with a synchronized
`temperature(nx,ny)` field. The normal predictor is selected independently from
the Riemann solver and supports PCM, frozen-composition characteristic PLM, or
time-traced characteristic PPM. For the y direction, momentum and primitive
velocity components are rotated into the x-normal ordering, evaluated by the
same predictor and selected HLLC/Rusanov/PeleC kernel, and rotated back.

```text
cell-centered conserved state
        ↓
NASA7 conserved-to-primitive recovery
        ↓
normal predictor selector
  ├─ PCM
  ├─ characteristic PLM
  └─ characteristic PPM
       ├─ five-point parabolic reconstruction
       ├─ u-c / u / u+c profile integration
       ├─ optional bounded contact steepening
       └─ optional shock flattening
        ↓
provisional selected x/y Riemann fluxes
        ↓
conservative transverse half-step correction
  U_face* = U_face - dt/(2 d_t) (F_t,hi - F_t,lo)
        ↓
EOS/positivity bisection on the complete conserved face state
        ↓
final selected directional Riemann fluxes
        ↓
unsplit two-dimensional conservative update
```

The transverse limiter acts on the complete conserved vector, including every
species density and total energy. Because each directional species-flux block
closes to the corresponding mass flux, the corrected face state retains
`sum(rho*Y_k)=rho` rather than repairing species independently after the hydro
correction. Both x-normal and y-normal characteristic-PLM/PPM reductions agree
with the corresponding 1D update at roundoff.

Chemistry uses the same cell-local adiabatic constant-volume solver as the 1D
path and is Strang split around the unsplit CTU hydro step. This path currently
requires periodic boundaries. The characteristic-PPM option supplies a
PeleC-style normal predictor in each coordinate direction before the existing
full-state CTU correction. Complete PeleC multidimensional PPM transverse/corner
tracing remains intentionally outside the current claim.

## Reaction-flow coupling

With molecular transport disabled, the coupled integrator retains the
reaction--hydro--reaction Strang sequence. With transport enabled, the symmetric
composition is:

```text
reaction(dt/2)
      ↓
transport(dt/2)
      ↓
hydro(dt)
      ↓
transport(dt/2)
      ↓
reaction(dt/2)
```

Each cell reaction solve holds density, all momentum components, and total-energy density fixed. Composition changes are written back to `rho*Y_k`, and temperature is recovered from the unchanged specific internal energy. This avoids adding a separate heat-release source on top of formation-energy-inclusive NASA7 internal energy.

## Verification separation

The architecture retains independent gates for:

1. constant-`gamma` hydro;
2. passive multispecies transport;
3. NASA7 thermodynamics;
4. zero-dimensional elementary chemistry;
5. composition-dependent reactive hydro;
6. reaction-flow splitting;
7. reactive two-dimensional CTU and dimensional reduction;
8. one- and two-dimensional molecular transport and Cantera coefficient
   qualification;
9. physical boundaries and full pressure-dependent H2/O2 chemistry;
10. MPI decomposition, multispecies hydro, implicit chemistry, molecular
    transport, and coupled reactive splitting;
11. dynamic reactive AMR hydro, chemistry, molecular transport, regridding,
    and conservative coarse/fine synchronization.

The homogeneous reactive field must reduce to independent zero-dimensional cell chemistry. The nonuniform hotspot must create finite pressure and velocity responses while preserving global mass, momentum, and total energy.

## Design constraints

1. Lower-order and constant-`gamma` baselines remain active in CI.
2. Invalid density, pressure, composition, temperature, or energy inversion fails explicitly.
3. Formation-energy offsets are retained in `h` and `u`.
4. Reverse rates use the same NASA7 records as the energy equation.
5. Species fluxes close exactly to the shared mass flux.
6. Chemistry does not independently modify `rhoE` in the adiabatic constant-volume substep.
7. The current characteristic basis assumes frozen composition across each acoustic solve.
8. Rusanov remains the robustness baseline. HLLC is the independent
   contact-resolving comparison; the PeleC-style path follows the upstream
   acoustic star-state and wave-interpolation sequence with the NASA7 EOS.
9. The four-reaction chemistry subset remains a lightweight regression path;
   the selectable ten-species, 29-reaction mechanism is the full H2/O2 path.
10. Contact steepening is explicitly bounded to half of the canonical detector
    strength until a complete general-EOS PPM/HLLC characteristic system is
    available.
11. The present transport layer excludes Soret, Dufour, multicomponent Stefan--
    Maxwell diffusion, polar corrections, and bulk viscosity.
12. Molecular transport is qualified in serial 1D/2D and distributed 1D paths;
    Soret and multicomponent diffusion remain excluded.

## Reactive PPM path

`reactive_1d_mod` keeps four independently selectable paths:

- `pcm`, the first-order robustness baseline;
- `characteristic_plm`, the frozen-composition MUSCL-Hancock path;
- `ppm`, the semidiscrete componentwise monotone path advanced by SSPRK3;
- `characteristic_ppm`, a time-centered normal predictor using PeleC's
  five-point parabolic reconstruction and `u-c`, `u`, `u+c` profile
  integration.

For multilevel AMR, `characteristic_ppm` can additionally enable the hybrid
WENO switch. `reconstruction_weno_mod` then replaces the PPM edge formula with
the selected WENO5-JS, WENO5-Z, WENO7-Z, or WENO3-Z nonlinear reconstruction;
the existing parabolic profile integration, frozen-composition characteristic
projection, SSPRK3 advancement, coarse-time ghost interpolation, and reflux
remain unchanged. Keeping WENO as an edge-reconstruction policy rather than a
fifth top-level hydro method follows the corresponding PeleC PPM control flow.

The characteristic PPM path carries species and transverse velocities on the
middle wave, projects density/normal velocity/pressure over the frozen mixture
acoustic basis, and converts the final face state through the NASA7 EOS. Its
one-dimensional shock-flattening coefficient follows PeleC `Godunov.H`.
Contact steepening is a separate Colella--Woodward-style detector applied to
density and species only. Both controls are opt-in and are rejected by the
configuration reader for other reconstruction modes.


## Reactive two-dimensional molecular transport

`reactive_transport_2d_mod` evaluates x/y face transport fluxes and advances
their conservative divergence independently of the CTU hydro operator.


## Physical boundary layer

`reactive_boundary_2d_mod` owns four typed faces and samples periodic, wall,
inflow, or outflow ghost states. Hydro and molecular transport use explicit
lower/upper face arrays. Solid walls receive a pressure-only inviscid flux,
mirrored velocity/temperature transport gradients, and zero species flux by
default. A wall face may instead own a prescribed species mass-flux vector.
The vector is oriented wall-to-gas, must sum to zero, and is translated to the
coordinate face orientation by `reactive_transport_2d_mod`; its species
enthalpy is added to total-energy flux and scaled with the same positivity
limiter as the species flux. This separates the boundary transport contract
from future catalytic surface kinetics.

## Full pressure-dependent H2/O2 chemistry

The reactive applications dispatch either the seven-species elementary model or a ten-species, 29-reaction model. The full path reuses the variable-width conserved state and advances each cell with the implicit constant-volume reactor.

## MPI one-dimensional path

`mpi_domain_1d_mod` owns uneven contiguous decomposition, periodic nonblocking
halo exchange, global reductions, and ordered `MPI_Gatherv` output. The MPI
drivers allocate only rank-local state plus two ghost cells.

```text
rank-local conserved state + temperature
        ↓
periodic state/temperature halo exchange
        ↓
global hydro/transport timestep reduction
        ↓
implicit chemistry(dt/2)
        ↓
SSPRK2 molecular transport(dt/2)
        ↓
conservative general-EOS Rusanov hydro(dt)
        ↓
SSPRK2 molecular transport(dt/2)
        ↓
implicit chemistry(dt/2)
```

Every coupled attempt is transactional. A failure on any rank is reduced across
the communicator, all ranks restore the pre-attempt state, and the scheduler
retries a smaller interval. CI compares complete gathered fields for 1, 2, 4,
and 8 ranks in both Debug and Release builds.

## AMR one-dimensional foundation

`amr_hierarchy_1d_mod` retains a reusable adjacent-level interface without
coupling AMR ownership to a particular fluid state width. Each child patch is
strictly nested inside its parent and uses an integer refinement ratio.

```text
coarse cell averages
        ↓ conservative limited prolongation
fine patch cell averages
        ↓ refinement-ratio subcycling
fine accumulated interface fluxes
        ↓ restriction + flux-register reflux
synchronized conservative composite state
```

Restriction replaces covered coarse cells with fine volume averages. Reflux
corrects the two uncovered coarse cells adjacent to the refined patch using the
time-integrated difference between fine and coarse interface fluxes.

The arbitrary-depth hierarchy owns an allocatable sequence of these adjacent
interfaces. Every interface may use a different refinement ratio, and each
level owns a separately allocated field because patch cell counts differ.
Initialization propagates physical bounds and spacing through the complete
nested chain. Composite integration counts each parent only outside its child,
then counts the deepest level in full. Synchronization applies reflux and
average-down from the deepest interface toward the root.

`amr_multilevel_reactive_1d_mod` adds separately allocated conserved-state and
temperature fields with one face-adjacent ghost and four PPM stencil layers at
every level. It recursively advances
each parent once and its child `r` times for hydro or `r^2` times for explicit
molecular transport, then synchronizes that relation before returning to the
next coarser caller. Chemistry is advanced on every level and averaged down
from deepest to root. The complete split update is transactional across the
hierarchy.

`amr_regrid_1d_mod` tags a selected state component using a normalized local
jump with an absolute floor, buffers the resulting tag interval, and constructs
one deterministic minimum-width patch. Regridding first averages old fine data
onto the coarse state, conservatively prolongs the new patch, then restores old
fine values on any same-resolution overlap. Consequently cells leaving a patch
retain their fine average, newly refined cells retain their coarse average, and
unchanged fine cells retain their full resolution. An empty tag set removes the
fine level after average-down. In the multilevel PPM/WENO engine, outflow
boundary tags may extend the patch to the matching physical edge. That side
uses fine-level constant-extrapolation ghosts; the opposite side retains
parent-interpolated coarse/fine ghosts. The legacy two-level PCM/PLM engine
continues to suppress boundary tags.

`amr_multipatch_1d_mod` represents an ordered set of separated child patches
over one parent level. It applies prolongation, average-down, reflux, and
composite integration across the set while excluding each covered parent
interval exactly once. `amr_regrid_1d_mod` can split disconnected tags into
candidate patches, expand each candidate, and coalesce candidates whose final
bounds touch or overlap. Regridding first synchronizes every old patch, builds
the new set, and restores all equal-resolution fine intersections even when a
patch moves, splits, or is repartitioned.

`amr_multipatch_reactive_1d_mod` qualifies fixed two-level reactive flow on
separated patches. The root hydro advances once, every child advances `r`
times from the same time-interpolated parent start/end states, and one register
per child is refluxed before set-wide average-down. Molecular transport uses
the same patch-wise synchronization with `r^2` child substeps per half
interval. Cell-local chemistry advances every root and child cell before a
set-wide average-down. The complete `R-T-H-T-R` interval is transactional and
reuses the existing characteristic PPM/WENO, diffusion, chemistry, and
coarse/fine ghost kernels.

The reactive driver retains the overlap-preserving two-level implementation for
PCM/PLM with `amr_max_levels = 2`. A larger level limit or either PPM option
selects the multilevel engine. It tags each parent, suppresses refinement at
interior patch boundaries, and constructs the next strictly nested child until
tags end or the configured depth is reached. PPM planning additionally reserves
the parent footprint required by four fine ghost layers and a limited parent
slope.
At a regrid point, it averages the old hierarchy deepest-to-root, plans a new
chain from the synchronized root, and rebuilds changed children by conservative
prolongation. It then maps every common old/new level in physical coordinates
and restores cells where spacing and cell boundaries align. A final
deepest-to-root average-down propagates retained fine information consistently.
Recursive output emits the left uncovered parent region, its
child, and the right uncovered region, producing ordered exact domain coverage.

The main dynamic reactive application selects the two-level patch-set engine
when `amr_multipatch_enabled` is true. It tags the synchronized root, clusters
disconnected features with a configurable maximum gap, constrains periodic
patches away from the domain seam, rebuilds the set at the configured interval,
retains every aligned old/new fine intersection, and emits ordered composite
CSV output.

`amr_patch_tree_1d_mod` composes those patch sets into an arbitrary-depth
forest. Each relation stores one child set per flattened parent patch and a
prefix-offset map from parent-local child indices to the next flattened level.
This permits a parent to own zero or more separated children while another
parent continues refining. Geometry validation reconstructs every parent
patch's physical extent, and field operations recursively prolong, average
down deepest-to-root, and integrate by replacing each covered parent interval
exactly once. Refinement ratios may differ between levels. A matching nested
register layout stores one flux register per child inside each parent-owned
set. Tree synchronization walks relations deepest-to-root and performs
set-wide reflux followed by average-down transactionally at every parent.

`amr_patch_tree_reactive_1d_mod` binds a reactive state, temperature, and wide
ghost storage to every flattened tree patch. A coarse interval is advanced by
a depth-first recursion: advance the parent once, accumulate its boundary
fluxes into every local child register, advance every child for `r` substeps
with time-interpolated parent ghosts, then reflux and average down that complete
child set. Independent branches may terminate at different levels. The public
advance is transactional and records the accepted advance count at every
level.

The public patch-tree step wraps that recursion in symmetric
chemistry--hydro--chemistry splitting. Chemistry advances every stored patch
for the same physical half interval, then average-down walks deepest-to-root so
covered parent cells again represent their children before the next operator.
The outer transaction restores state, temperature, time, and counters together
if either reaction half-step or hydro fails.

Molecular transport uses the same parent-owned register topology with
parabolic subcycling. A parent advances once for its transport interval; each
child advances `r^2` times with midpoint-interpolated parent ghosts, returns its
time-integrated diffusive boundary fluxes, and is refluxed and averaged down
before returning. The timestep reduction applies the square of every
cumulative refinement ratio to each fine-patch transport limit.

Runtime regridding first synchronizes the old tree to the root, constructs a
new tree from an explicit branching plan, and conservatively prolongs that
root. For every common level with matching spacing, all old/new patch pairs are
intersected in physical coordinates and aligned cells are copied exactly even
when their parent ownership changes. A final deepest-to-root average-down
restores covered-parent consistency. State, temperature, time, advance
counters, and regrid statistics roll back together on failure.

Automatic planning uses that synchronized root as its deterministic source.
At each prospective relation it tags every parent patch independently,
restricts tags to cells that retain enough parent stencil support, and
clusters disconnected features with the configured gap, buffer, and minimum
width. Parent-local collections flatten in parent order into the next level.
The resulting children are prolonged from the root and become the parents for
the next tagging pass until tags end or `amr_max_levels` is reached. The plan
then enters the same transactional rebuild and overlap-transfer path as an
explicit plan.

Patch-tree child sets may retain independently owned adjacent intervals. At
each child substep the parent fills coarse/fine ghosts first, then sibling
interiors overwrite every covered face ghost and available PPM/WENO wide
layer by global fine index. After the siblings advance, their two returned
time-integrated fluxes at each shared face are replaced by one arithmetic-mean
flux. Conservative corrections apply that owned flux to the two adjacent fine
cells. The internal face is suppressed in both sibling flux registers, so
only genuine coarse/fine sides participate in reflux. The same sequence is
used for hyperbolic and molecular-transport recursion.

The qualified patch-tree path covers interior separated or adjacent patches,
PCM/PPM hydro, elementary chemistry, molecular transport, and explicit or
tag-driven runtime rebuilds. One-sided periodic-seam refinement remains
separate.

`mpi_amr_patch_1d_mod` adds the first distribution boundary without importing
MPI into the serial AMR modules. Every rank first proves that the integer patch
topology and physical root extent are identical. A deterministic greedy work
schedule then gives each patch exactly one owner; ties go to the lowest rank.
For level `l`, patch work is its cell count times the cumulative refinement
ratio raised to exponent 0, 1, or 2. These select storage/cell weighting,
hyperbolic `r` subcycling, or parabolic `r^2` subcycling respectively. The
work model and 64-bit per-patch/per-rank totals are validated collectively and
preserved through explicit and tag-driven sparse regrids. Generic patch fields
remain allocated on every rank in this
bridge, but only the owner is authoritative. Collective broadcasts refresh
the replicas, while adjacent sibling faces broadcast the owner's boundary
cells into explicit left/right halo objects for up to four stencil layers.

For sparse fields, `sparse_patch_tree_reactive_timestep_1d` evaluates
hyperbolic and optional parabolic stability limits only on allocated owner
payloads. Each level-local limit is converted to its root-step equivalent with
cumulative `r` or `r^2` scaling before one communicator-wide minimum
reduction. Ranks with no local patches contribute the neutral `huge()` value,
while any invalid owner state or missing transport database is rejected
collectively before a timestep is published.

This separation establishes rank ownership and communication ordering before
changing the recursive reactive integrator. In `0.50.0`, the chemistry
operator uses that ownership directly: every rank walks the same patch order,
only the owner integrates a patch, all ranks reduce the local success flag,
and the owner broadcasts the accepted complete reactive patch. Once all
patches finish, every replica performs the same deepest-to-root average-down,
temperature recovery, and ghost refresh. A backup taken after initial owner
synchronization makes any rejected owner update a communicator-wide rollback.

In `0.51.0`, serial and MPI recursion share a one-patch hydro kernel. Every
rank traverses the same parent, child, and substep order, but only the patch
owner executes that kernel. Collective acceptance precedes broadcasts of the
owner's interval-start state, complete face-flux field, accepted patch state,
and level counter. Synchronized replicas then apply the existing
time-interpolated child ghost fill, adjacent-sibling exchange and
time-integrated flux reconciliation, coarse/fine flux-register accumulation,
reflux, average-down, temperature recovery, and final ghost refresh. The
owner-synchronized backup makes rejection at any recursion depth a global
transactional rollback.

This is an owner-authoritative replicated bridge, not a sparse distributed
stage implementation: parent/child and fine/fine operations still execute on
every replica after owner broadcasts. Owner-only molecular transport, sparse
rank-local storage, migration after regrid, stage-synchronous point-to-point
halos, and scalable communication schedules remain outside the `0.51.0`
qualification boundary.

In `0.52.0`, molecular transport gains the same shared-kernel boundary and
owner recursion. The owner executes the complete SSPRK2 viscous, conductive,
and mixture-averaged diffusion update and broadcasts its interval-start state,
effective face fluxes, accepted patch, and transport counter. Child recursion
uses the serial parabolic schedule of `r^2` substeps per relation, including
time-interpolated parent ghosts, adjacent shared diffusive fluxes,
coarse/fine registers, reflux, average-down, and temperature recovery.
Collective rejection restores the pre-transport replica on every rank.

Chemistry, hydro, and transport now each have qualified owner-only entry
points. A single outer transaction composing all three operators, sparse
rank-local patch allocation, migration after regrid, and scalable
point-to-point communication remain outside the `0.52.0` boundary.

In `0.53.0`, `advance_owned_patch_tree_reactive_1d` composes the three owner
operators as `R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)`. Before taking the outer
backup, owner synchronization now includes the complete patch fields plus
time, step, hydro and transport counters, and regrid statistics, with the root
patch owner authoritative for global bookkeeping. Every stage retains its
inner collective acceptance, while the outer wrapper restores the common
pre-step tree and reports zero accepted calls if any later stage fails.

This establishes one qualified distributed full-physics transaction with the
same operator order and final ghost contract as serial patch-tree AMR. Patch
arrays and hierarchy operations remain replicated; sparse rank-local
allocation, migration after regrid, and scalable point-to-point communication
remain outside the `0.53.0` boundary.

In `0.54.0`, `mpi_amr_sparse_patch_1d_mod` introduces a rank-local reactive
patch container. Hierarchy and owner-map descriptors remain replicated, but
the six allocatable field payloads are present only for patches owned by the
local rank. Scatter first resolves the owner-authoritative replicated tree;
gather writes local owners into a supplied replica and reuses the qualified
owner synchronization to reconstruct every patch and global bookkeeping.

Same-hierarchy owner-map changes migrate one complete patch at a time from
the old owner. The current transition uses a collective broadcast as a
correctness-first bridge, although only the new owner retains the payload.
Physics still executes through the replicated `0.53.0` entry points. Direct
sparse physics, topology-changing regrid transfer, and point-to-point
communication schedules remain outside the `0.54.0` boundary.

In `0.55.0`, `advance_sparse_patch_tree_chemistry_1d` advances only locally
allocated patch payloads. Deepest-to-root synchronization streams one child
interior at a time to the parent owner for conservative average-down and
temperature recovery. Root physical ghosts and parent/fine ghosts are updated
only on their owners; adjacent sibling state is likewise streamed one source
patch at a time so normal and PPM-wide ghosts can be replaced locally.

Collective acceptance surrounds every patch and synchronization boundary. A
failure restores each rank's sparse backup and reports zero accepted calls.
Hydro, molecular transport, and the combined `R-T-H-T-R` transaction still use
the replicated bridge. Broadcast-based synchronization remains a temporary
correctness schedule outside the final point-to-point design.

In `0.56.0`, `advance_sparse_patch_tree_hydro_1d` runs the finite-volume patch
kernel only on the sparse owner. Each recursion broadcasts that patch's
interval-start state, accepted interval-end state, and effective face flux so
child owners can apply time-interpolated ghost data and every rank can update
the same small flux-register metadata. Child recursion, mixed-ratio
subcycling, and adjacent-face integral reconciliation follow the serial order.

After child subcycles, interiors are streamed to the parent owner for reflux,
average-down, and temperature recovery. Only the owner mutates the persistent
parent payload. A sparse final ghost refresh and global time/step update close
the transaction; rejection restores rank-local payloads and all counters.
Molecular transport and combined full physics remain replicated, and the
communication schedule is still broadcast-based.

In `0.57.0`, `advance_sparse_patch_tree_transport_1d` applies the same sparse
recursion boundary to SSPRK2 molecular transport. The owner computes each
viscous, conductive, and mixture-diffusion patch interval. Parent interval
states and effective diffusive fluxes are streamed for `r²` child subcycles,
and adjacent time-integrated diffusive faces are reconciled on their owners.

Compact diffusive registers remain replicated; child interiors return to the
parent owner for reflux, average-down, and temperature recovery. Final ghosts
and rollback remain sparse. All three component operators now have direct
sparse entry points, but their combined `R-T-H-T-R` transaction and the
communication schedule remain future work.

In `0.58.0`, `advance_sparse_patch_tree_reactive_1d` composes those direct
sparse operators as `R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)`. It validates the
optional transport database collectively before mutation, takes one outer
sparse backup, and accepts every stage across the communicator. A rejected
later stage restores fields, ghosts, time, steps, and both level-counter
vectors on every rank and reports zero committed operator calls.

The normal distributed full-physics path can therefore stay rank-local from
entry through final ghost refresh and validity checking. Each component still
uses its own transactional backup and correctness-first collective streaming;
topology-changing regrid transfer and point-to-point schedules remain outside
the `0.58.0` qualification boundary.

In `0.59.0`, `regrid_sparse_patch_tree_reactive_1d` adds an explicit-plan
topology transition around the serial qualified rebuild. It materializes the
old sparse tree collectively, applies conservative average-down, root
prolongation, and same-resolution overlap transfer, verifies identical change
decisions and transfer counts across ranks, then constructs a new deterministic
owner map and scatters only the rebuilt owners.

An unchanged plan increments only the regrid-evaluation counter and retains
the existing distribution. A changed plan commits the rebuilt hierarchy,
payloads, counters, and distribution together; any invalid plan or collective
failure restores the original sparse solution and owner map. Persistent state
is again globally single-copy after return, although this correctness-first
transition temporarily materializes a replica. Direct tag planning and
point-to-point regrid transfer remain outside the `0.59.0` boundary.

In `0.60.0`, `regrid_tagged_sparse_patch_tree_reactive_1d` materializes the
current owner state once, applies parent-local gradient tagging and clustering
through the configured maximum depth, and sends the resulting serial-qualified
tree through the shared sparse regrid commit helper. Tagged-cell, change, and
overlap counts must agree across every rank before the solution and owner map
can commit.

Disconnected root features can therefore generate separate children through
four levels from the sparse public API. Repeating the same tag decision is a
no-op except for evaluation accounting, and an invalid tag component restores
the original sparse solution and distribution. Tag planning still operates on
the temporary correctness replica; owner-local tag construction and
point-to-point topology transfer remain outside the `0.60.0` boundary.

In `0.61.0`, same-hierarchy ownership migration no longer broadcasts every
patch to every rank. For each changed owner, the old owner packs state,
temperature, narrow ghosts, and wide ghosts into one contiguous payload and
sends it directly to the new owner. Unchanged owners copy locally, and ranks
unrelated to that patch carry no payload buffer.

A collective acceptance boundary follows each ordered patch transfer, so the
new sparse container is published only after all direct messages and final
shape checks succeed. This removes replicated migration traffic while keeping
the deterministic correctness schedule. Physics-stage streaming and topology
regrid remain collective outside the `0.61.0` boundary.

In `0.62.0`, adjacent sparse siblings no longer broadcast each complete source
patch through the communicator. Each adjacent pair packs only the state and
temperature boundary layers required by the active reconstruction: one layer
for narrow ghosts or four layers for PPM. Different owners exchange the two
directional payloads with one blocking `MPI_Sendrecv`; the left owner counts
the completed pair once, while ranks unrelated to the face allocate nothing.

Siblings on the same owner copy their boundary layers locally. Both paths
refresh the narrow state/temperature ghost and the optional wide PPM arrays,
preserving the existing recursive chemistry, hydro, and transport results.
Parent interval streaming, flux reconciliation, average-down, and
topology-changing regrid transfer remain collective outside the `0.62.0`
boundary.

In `0.63.0`, child interiors consumed by chemistry average-down and by
hydro/transport reflux synchronization no longer broadcast from the child
owner. A shared transfer helper sends the contiguous interior state directly
to the corresponding parent owner. When both patches have the same owner, the
state is copied locally; unrelated ranks allocate no child payload.

The deterministic child order keeps blocking sends and receives matched across
the recursive traversal. Chemistry reports completed remote transfers on child
owners, and the communicator sum must equal an independently derived
cross-owner child count. Parent interval and flux streaming, parent-to-child
ghost fill, flux reconciliation, and topology-changing regrid transfer remain
collective outside the `0.63.0` boundary.

In `0.64.0`, final sparse ghost refresh no longer broadcasts a complete parent
state to the communicator. The parent owner derives the distinct ranks that
own at least one of its children and sends the parent state once to each remote
recipient. A child owner receives one reusable copy even when it owns several
children of that parent, and ranks without a related parent or child allocate
no state buffer.

Same-owner child fill reads the local parent directly. The parent owner counts
each successful remote-recipient send, and the communicator sum must equal the
owner-map-derived number of distinct remote child owners. Recursive hydro and
transport interval states, face fluxes, bookkeeping counters, flux
reconciliation, and topology-changing regrid transfer remain collective
outside the `0.64.0` boundary.

In `0.65.0`, recursive sparse hydro and molecular transport remove their final
owner broadcasts. Each patch owner keeps the full update flux and its flux
register authoritative locally. It packs interval start/end states once for
each distinct remote child owner, while a child patch sends only its two
time-integrated boundary flux vectors to its parent owner.

The parent owner computes adjacent shared-face fluxes and sends each resulting
cell correction only when the affected child has a different owner. Coarse and
fine register accumulation and final reflux therefore stay on the parent
owner. Hydro and transport level counters accumulate owner-local deltas during
recursion and use one array reduction at the stage boundary, replacing the
per-patch counter broadcasts.

The sparse physics module now contains no `MPI_Bcast`. Exact global counts for
interval-state fanout, child boundary-flux returns, and shared-flux corrections
are derived independently from hierarchy ownership and subcycle weights for
both hyperbolic `r` and parabolic `r^2` recursion. Topology-changing regrid
transfer and its temporary correctness replica remain outside the `0.65.0`
boundary.

In `0.66.0`, the explicit-plan sparse regrid path no longer materializes either
the old or rebuilt field tree on every rank. All ranks construct only the
compact candidate hierarchy and deterministic owner map. Starting at the root,
each parent owner applies the qualified conservative prolongation and sends one
interior-state payload only to a remote child owner; same-owner children remain
local.

After owner-local child allocation, every same-resolution old/new patch overlap
is enumerated from replicated geometry metadata. State and temperature cells
move directly from the old patch owner to the new patch owner, while local
overlaps copy without MPI. Deepest-to-root sparse average-down and the existing
point-to-point ghost refresh complete the new hierarchy before transactional
publication.

The application independently derives exact global counts for cross-owner
prolongation and overlap messages and retains exact serial full-field parity.
An unchanged explicit plan only increments evaluation accounting. Tag-driven
plan construction still materializes a temporary correctness replica and is
outside the `0.66.0` boundary.

In `0.67.0`, tag-driven regrid also remains field-sparse. A rank-local planning
copy first performs the qualified distributed average-down. At each candidate
depth, only a parent owner evaluates its conserved-state gradients and clusters
tagged cells. Two integer reductions publish tagged counts and child bounds as
compact topology metadata; no state or temperature field is replicated.

Candidate child states are conservatively prolongated through the same direct
owner path as explicit regrid. Exact counts cover parent-owner tag evaluations,
candidate prolongation messages, final topology prolongation, and retained
overlap transfers. The final plan enters the `0.66.0` direct transactional
rebuild, so both public topology-change APIs now preserve globally single-copy
field storage throughout.

Installing a new conserved field resets the temperature recovery seed before
EOS inversion. This makes rebuilt temperatures independent of stale guesses and
preserves exact serial/distributed parity even when a tagged conserved state was
edited without updating its cached temperature.

The `0.73.0` public driver adds a rank-neutral persistence boundary. At a
configured coarse-step cadence, all owners participate in the existing
materialization operation, rank zero writes a versioned patch-tree checkpoint,
and the sparse solution remains authoritative if the run continues. Restart
loads the hierarchy and fields without any owner identifiers, validates the
mechanism layout and root geometry, and schedules the recovered patches over
the current communicator before normal sparse advancement resumes. Thus
changing the MPI rank count changes placement, not the persisted numerical
state.

## Reactive AMR time advancement

`amr_reactive_1d_mod` owns a coarse reactive state, an optional fine state, both
temperature fields, hierarchy metadata, simulation time, and regrid counters.
One accepted coarse interval with all operators enabled is:

```text
chemistry(dt/2) on coarse and fine
        ↓
SSPRK2 molecular transport(dt/2) on coarse and fine
  with r^2 fine diffusion substeps, reflux, and average-down
        ↓
coarse PCM/PLM/PPM hydro(dt) and interface-flux capture
        ↓
fine PCM/PLM/PPM hydro(dt/r), repeated r times
  with time-interpolated coarse ghost states
        ↓
flux-register reflux + average-down
        ↓
SSPRK2 molecular transport(dt/2) on coarse and fine
  with r^2 fine diffusion substeps, reflux, and average-down
        ↓
chemistry(dt/2) on coarse and fine + average-down
        ↓
optional tagging and conservative regrid
```

The coarse hyperbolic CFL limit is combined with `r` times the fine hyperbolic
limit. The parabolic limit is combined with `r^2` times the fine transport
limit, so every advective and diffusive substep is stable. A complete solution
copy makes the interval
transactional: any EOS, Riemann, chemistry, transfer, or synchronization
failure restores both levels and all hierarchy metadata. Composite output emits
uncovered coarse cells and fine cells exactly once and in coordinate order.

The optional AMR PLM path reconstructs density, velocity, pressure, and mass
fractions with the configured MC or minmod limiter. Face mass fractions are
clipped and renormalized before conversion back to the conserved general-EOS
state. SSPRK2 evaluates two flux divergences and returns their arithmetic mean
to the flux register, so the reflux correction represents the same conservative
update applied to each level. Fine PLM substeps hold a coarse-time midpoint
ghost state during both SSPRK2 stages.

The AMR PPM path reconstructs both cells adjacent to every coarse/fine face.
Four exterior conserved states and temperatures cover the widest
characteristic stencil. Fine ghost cell centers are obtained by linearly
interpolating parent conserved averages in time and applying an MC-limited
conservative parent slope in space. Each fine substep holds its midpoint ghost
data through an SSPRK3 update. The returned interface flux is
`(F0 + F1)/6 + 2 F2/3`, exactly matching the conservative SSPRK3 state update,
so reflux uses the flux that advanced the level.

The AMR transport path evaluates the same Newtonian stress, Fourier heat flux,
mixture-averaged species diffusion, barodiffusion, correction velocity, and
species-enthalpy flux as the uniform 1D path. Coarse and fine levels each use
SSPRK2, and the mean stage flux is accumulated in a separate diffusive flux
register. Fine transport uses `r^2` substeps over each half interval. At a
coarse/fine face, gradients use the actual distance between the adjacent
coarse and fine cell centers. Reflux and average-down complete each transport
half step before the next split operator begins.

## Embedded-boundary geometry foundation

The `0.74.0` EB subsystem is independent of the regular-cell hydro path. A
caller supplies finite level-set values on Cartesian nodes, with positive
values denoting fluid. `eb_geometry_2d_mod` allocates cell-volume fractions,
x-face and y-face open-area fractions, and regular/cut/covered classifications.

Each quadrilateral follows one fixed lower-left to upper-right diagonal. The
level set is affine on both resulting triangles, and each positive polygon is
clipped and integrated in normalized cell coordinates. Shared faces use the
same endpoint interpolation once, so neighboring cells reference one open
fraction. Interface length, centroid, representative unit normal, and
integrated normal support a pressure-only reactive slip wall. The EB operator
combines open-face and wall fluxes, then the first-order redistribution module
stabilizes its right-hand side and transactionally recovers the updated
reactive state. A second path applies weighted StateRedist to the provisional
conserved state. It derives up to three merge neighbors from the
face-aperture normal, represents overlapping neighborhoods explicitly, and
uses partitioned self/neighbor weights before the EOS-validated commit. The
geometry also stores normalized fluid-volume centroids. With
`state_redist_max_order=2`, `Qhat` is located at the matching weighted
neighborhood centroid, a connected 3-by-3 or active grown 5-by-5 least-squares
fit supplies linear slopes, and centroid plus merge-recipient limiters bound
the reconstructed values without changing their zero first moment. Flux
construction is connected by `eb_reactive_hydro_2d_mod` and
`eb_reactive_reconstruction_2d_mod`. PCM remains the baseline; selectable
characteristic PLM forms frozen-composition limited slopes only where both
normal neighbors are active, traces the normal waves, and falls back to zero
slope beside covered cells and outer boundaries. Riemann fluxes are first
formed at Cartesian face centers and then linearly interpolated in the
tangential direction to normalized open-face centroids. Zero-gradient domain
faces reuse their adjacent fluid cell, and the divergence feeds the weighted
StateRedist transaction. `reactive_eb_2d_driver_mod` turns this
operator into a standalone time-dependent application: it builds a plane or
circle level set from the input file, evaluates the general-EOS CFL rate over
active cells only, advances to the clipped final time, and emits
volume-fraction-weighted diagnostics and geometry-aware CSV output. The public
driver accepts its qualified PCM-or-PLM/outflow contract and can wrap hydro in
active-cell reaction half steps. Covered cells are excluded from both reactor
calls, and candidate arrays make the complete reaction--hydro--reaction step
transactional. Molecular transport remains rejected.

`amr_eb_hierarchy_2d_mod` provides the first bridge from this EB geometry to
the AMR data model. It describes one aligned rectangular level-two patch by its
parent bounds and integer refinement ratio. Construction rejects inconsistent
physical bounds, spacing, dimensions, or parent/child volume measures. The
generic restriction replaces only covered parent-patch cells with fine
fluid-volume-weighted states; a composite integral omits that parent region and
counts the fine patch instead. Its reactive wrapper treats state plus recovered
temperature as one transaction and retains covered-parent data.

`amr_eb_flux_register_2d_mod` adds the matching conservative interface path.
The register stores accumulated state corrections on coarse cells immediately
outside the patch. Coarse faces contribute their open-area flux with the
outward finite-volume sign; fine faces contribute the opposite time-integrated
sum over aligned subfaces. Reflux applies regular-cell corrections directly.
For a cut interface cell it follows AMReX EB re-reflux: keep the `kappa` share,
scatter the remaining share over connected neighbors by their fluid volumes,
and transfer recipients covered by the fine rectangle to the corresponding
fine children. Candidate registers and level arrays make reactive state and
temperature recovery atomic.

`amr_eb_reactive_2d_mod` composes those kernels into one static two-level
hydrodynamic interval. Piecewise-constant prolongation injects each parent
state into its `r` by `r` children and recovers active-child temperature. The
coarse level advances once, and the fine patch advances `r` times with
`dt_f=dt_c/r`. At every fine substep, its four exterior face-state arrays are
filled from adjacent coarse cells at an interpolated coarse time. The exact
centroid fluxes used by both level advances accumulate into the register;
reactive re-reflux and EB average-down complete a single all-or-nothing state
and temperature update.

`simulation_config_reactive_eb_amr_2d_mod`,
`reactive_eb_amr_2d_driver_mod`, and `pelef_reactive_eb_amr_2d` provide the
first runnable hierarchy shell. A third namelist group supplies inclusive
coarse patch bounds, refinement ratio, fine output path, and optional dynamic
regridding controls. Initialization builds the same configured level set
independently on both resolutions, checks the AMR geometry contract, and
prolongs the reactive state. Every coarse timestep is the smaller of the root
stability limit and `r` times the fine stability limit, then clipped to the
requested final time. The app writes the synchronized parent and child fields
separately after the time loop.

`amr_eb_regrid_2d_mod` adds the single-patch topology transaction. Four-neighbor
temperature jumps tag active internal root cells through combined relative and
absolute thresholds. A buffered tag bounding box is clamped to the strictly
internal region and grown to configured minimum extents. Before replacing the
patch, the old fine data are volume-weighted into the root. PCM initializes the
new fine rectangle, then matching global fine indices restore exact old overlap.
The new hierarchy is committed only after all active new fine cells pass EOS
recovery. The driver invokes this planner optionally at initialization and at a
configured accepted-step cadence. An unchanged plan retains the current patch.
With fine-patch removal enabled, an empty plan average-downs the complete child,
releases its state, temperature, geometry, and patch metadata, and leaves one
root level. A later active plan rebuilds the fine geometry and initializes it
by PCM from the synchronized root. The timestep and advance dispatch select the
two-level or single-level EB path from this lifecycle state, and inactive fine
output is omitted.

In `0.91.0`, `reactive_eb_amr_2d_driver_mod` composes chemistry around that
lifecycle-aware hydro dispatch. When a fine patch is active, candidate coarse
and fine states each receive a masked reaction half-step, followed by the
existing subcycled EB hydro/reflux transaction, a second reaction half-step on
both levels, and reactive average-down. The entire hierarchy remains private
until all stages succeed. With no fine patch, the same driver calls the
qualified single-level EB Strang operator. The application loads either the
elementary or full H2/O2 reaction set through the existing mechanism path.

In `0.92.0`, the same serial driver owns a versioned formatted checkpoint
schema. It stores conserved state and temperature on the root and, when active,
the fine rectangle; actual patch bounds; lifecycle/run counters; species names;
and a compatibility signature for mesh, EB geometry, chemistry, hydro,
redistribution, and regrid settings. Restart reconstructs geometry from the
current input rather than trusting serialized metrics, reads into private
candidates, recovers every active temperature from conserved state, and commits
only after the end marker. A root-only checkpoint keeps fine arrays and metadata
unallocated. Final time, step budget, output paths, and checkpoint scheduling
are intentionally restart-mutable.

In `0.93.0`, `amr_eb_regrid_2d_mod` adds a separate two-level multipatch
kernel. A deterministic flood fill clusters disconnected root tags with a
configurable maximum gap. Each cluster is buffered and expanded independently;
candidate rectangles whose two-cell separation would violate the EB
redistribution neighborhood are coalesced. The resulting collection owns an
ordered set of fine geometries, patch metadata, conserved states, and
temperatures. Topology replacement first averages the old set into a private
root, builds every new child by PCM, and then copies exact same-resolution data
over every old/new patch intersection before an EOS-validated commit.

The multipatch hydro transaction advances the root exactly once and advances
each child through `r` substeps using coarse-time-interpolated exterior states.
Each child owns a distinct EB flux register. The transaction refluxes the
children into a private root in order and then average-downs the full set. A
matching driver transaction applies masked reaction half-steps to the root and
all children around that hydro update, synchronizes after the second reaction,
and exposes no partial hierarchy on failure.

In `0.94.0`, the public EB AMR executable dispatches to this representation
when `multipatch_enabled` is true. A configured rectangle seeds the initial
set, after which the collection planner may replace it at initialization and
at the accepted-step regrid cadence. The timestep routine reduces the active
root CFL limit with `r` times every child limit. Each accepted interval calls
the set-wide Strang transaction, and a committed periodic regrid replaces the
root and complete child collection together. Output iterates children in
planner order and adds a stable patch number before the CSV extension.

In `0.94.0`, the version-one formatted checkpoint had one optional child
payload, so configuration and runtime validation rejected checkpoint or
restart controls whenever multipatch mode was enabled.

`0.95.0` adds a distinct patch-set checkpoint magic and schema without changing
that established single-patch format. The patch-set manifest stores the child
count and then each deterministic child's actual coarse bounds, dimensions,
conserved state, and temperature after the root and compatibility signature.
The reader rebuilds all geometries and states privately, recovers active
temperatures from conserved variables through the EOS, validates set-wide
separation and the terminal marker, and publishes the root and complete child
collection together. Scheduled writes occur after physics and periodic regrid
transactions, so restart resumes the accepted-step cadence.

In `0.96.0`, the configured single-patch path permits a rectangle to coincide
with a root physical side. Geometry construction uses the exact domain edge.
During every fine substep, a coarse/fine side still receives the established
coarse-time-interpolated exterior state, while a coincident outflow side uses
the current fine boundary cell as its zero-gradient exterior state. The flux
register already omits any side without an uncovered coarse/fine interface, so
no correction is accumulated or applied at the physical boundary. Fine state
and temperature are passed explicitly to exterior construction and validated
before the substep publishes any candidate hierarchy.

In `0.97.0`, temperature-gradient tagging visits every active root cell. At a
physical side it compares only in-domain active neighbors, giving a one-sided
gradient without fabricating a root ghost value. Single- and multipatch plan
bounds, buffer expansion, minimum-size growth, and component flood fill clamp
to indices `1:nx` and `1:ny`. The patch-set hydro transaction accepts those
domain-inclusive children and reuses the 0.96 physical-side exterior closure;
its flux registers continue to act only on true coarse/fine interfaces.
Topology replacement, child ordering, rollback, and collection separation are
unchanged.

In `0.98.0`, `amr_eb_multilevel_2d_mod` composes two aligned EB patch
descriptors into a strictly nested root/middle/finest hierarchy. Composite
integration owns root cells outside the middle patch, middle cells outside the
finest patch, and every finest cell exactly once. Synchronization first
average-downs finest into a private middle candidate, then that candidate into
a private root candidate. The reactive path recovers EOS-consistent
temperatures after both restrictions and publishes neither parent if either
stage fails.

In `0.99.0`, `amr_eb_multilevel_reactive_2d_mod` recursively advances that
hierarchy. The root advances once over `dt`; the middle advances `r1` times;
inside each middle interval the finest advances `r2` times. Each parent/child
pair owns an independent EB flux register. Finest reflux and average-down
finish every middle interval, then the accumulated middle flux refluxes the
root and a final deepest-first synchronization commits all three levels. The
finest patch must remain two middle cells from the middle boundary and its
coarse/fine interface must be fully regular. These checks reject the known
unsupported EB-cut nested-interface case before any state is published.

In `0.100.0`, the regular-interface restriction is replaced by a conservative
EB-cut closure. Before each middle update, the module records the authoritative
middle/finest composite integral. After inner reflux and average-down it adds
the time-integrated flux through the middle domain boundary to form the
expected integral. Any remaining density, total-energy, and species residual
is spread uniformly per fluid volume over active middle cells outside the
finest patch. Species corrections close to the density correction, and every
recipient temperature is recovered through the EOS before publication.

In `0.101.0`, `reactive_eb_amr_2d_driver_mod` composes the qualified
active-cell reactor with the static three-level hydro transaction. Reaction
half-steps advance private root, middle, and finest candidates before and
after recursive hydro. The second half-step is followed by reactive
finest-to-middle-to-root average-down, so chemistry applied independently on
overlapped parents cannot become the final authoritative state. Failure in
chemistry, hydro, EB-cut conservation closure, or EOS synchronization returns
all three original state and temperature fields.

In `0.102.0`, the public serial application can select a static three-level
mode. The root rectangle in the existing `eb_amr` namelist defines the middle
mesh; a second rectangle in middle indices defines the finest mesh and must
retain the qualified two-cell margin. Initialization prolongs root to middle
and middle to finest. Every accepted root interval selects
`min(dt0, r1*dt1, r1*r2*dt2)`, clips it to the final time, and commits one
three-level Strang transaction. Root, middle, and finest CSV files are written
only after successful completion.

In `0.103.0`, that lifecycle owns a separate three-level checkpoint schema.
The writer records the complete root, middle, and finest state and temperature
fields, the two nested rectangles and refinement ratios, accepted time and
step metadata, ordered mechanism layout, and the numerical compatibility
signature. The reader rebuilds all three EB geometries and validates the
entire stream, including its terminal marker, in private candidates. It then
recovers every active temperature through the EOS before publishing the
restored hierarchy. Scheduled and final writes occur only after accepted
three-level transactions; stop-after-write and restart preserve timestep
cadence without modifying the older checkpoint schemas.

In `0.104.0`, `dynamic_regridding` may instead keep the configured middle
level fixed and rebuild the finest rectangle from middle-level temperature
gradients. The planner operates only on the interior obtained by removing two
middle cells from every side. A topology transaction first average-downs the
old finest state, constructs and prolongs the tagged replacement, restores
overlapping fine cells exactly, validates EOS state, and then publishes the
new middle/finest pair. Initial and cadence-triggered regrids use the same
operation, and the public lifecycle reports committed topology changes.

In `0.105.0`, dynamic three-level mode selects a separate formatted checkpoint
magic and schema. The stream records the committed middle-to-finest bounds,
regrid count, interval, tagging thresholds, buffer and minimum-size controls
in addition to the established mechanism, physics, EB, state and time data.
Restart rebuilds the fixed root/middle hierarchy from configuration and the
finest geometry from the stored bounds, recovers all temperatures through the
EOS, and publishes only after the end marker and every compatibility check
succeed. Restored step and regrid counts preserve accepted-step cadence.

In `0.106.0`, the single-level reactive EB driver accepts the established
mixture molecular-transport database. Regular Cartesian diffusive fluxes are
evaluated with the configured physical boundaries and interpolated to open
face centroids. An EB-aware limiter budgets outgoing species against each
cell's fluid-volume inventory and scales the complete coupled face flux. The
open-area divergence has no embedded-face contribution, which defines an
adiabatic slip and species-impermeable wall for this milestone. Each SSPRK2
transport stage uses the existing StateRedist/EOS transaction, and the full
step composes `R/2 -> T/2 -> H -> T/2 -> R/2` without publishing partial work.

In `0.107.0`, `amr_eb_transport_2d_mod` lifts that operator into the
single-patch two-level hierarchy. A coarse transport Euler stage exposes its
EB face-centroid fluxes, while each ratio-subcycled fine stage samples
time-interpolated coarse exterior states. The existing EB flux register
accumulates both levels' time-integrated diffusive fluxes, then reactive reflux
and average-down synchronize the hierarchy. Two complete synchronized Euler
transactions form SSPRK2, so each stage is conservative before the RK average.
The public driver combines root and ratio-scaled fine transport limits and
composes this hierarchy operator as `R/2 -> T/2 -> H -> T/2 -> R/2`.

In `0.108.0`, `amr_eb_multilevel_transport_2d_mod` applies the same transport
transaction recursively to a root, middle, and finest hierarchy. Each middle
substep owns a complete finest subcycle and closes its inner flux register
before contributing middle flux to the outer register. EB-cut finest
interfaces reuse the conservative residual closure, and each synchronized
Euler stage restricts finest to middle before middle to root. The public
driver selects the minimum root-equivalent stability limit across all three
levels and retains the symmetric `R-T-H-T-R` composition.

Unsplit transverse prediction, fourth-order StateRedist slopes, periodic ghost
neighborhoods, thermal/viscous/catalytic walls, coarse-to-fine spatial slopes,
multipatch EB AMR molecular transport, locally resolved
PeleC-style multilevel redistribution, arbitrary depth, dynamic root/middle lifecycle ownership,
non-outflow refined boundaries, and EB AMR/MPI ownership remain outside this
subsystem. Dynamic three-level mode changes only the finest patch inside a
fixed middle level and rejects finest removal and siblings.
