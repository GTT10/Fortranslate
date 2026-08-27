# PeleF architecture

## Executable split

PeleF exposes twelve serial verification drivers and eight optional MPI drivers
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

In `0.109.0`, `amr_eb_multipatch_transport_2d_mod` advances the root
transport stage once and gives every separated sibling patch an independent
ratio-subcycled fine transaction and diffusive flux register. All children
sample the same time-interpolated root candidate. Their disjoint reflux
updates are accumulated before one patch-set average-down. If an embedded
boundary crosses any child interface, one set-wide composite residual closure
corrects active unrefined root recipients, closes density against species, and
recovers every corrected temperature through the EOS. Two complete
synchronized Euler transactions form SSPRK2, and the public driver includes
all child parabolic limits in the root timestep.

In `0.110.0`, `mpi_amr_eb_patch_2d_mod` introduces a distribution boundary
without importing MPI into the serial EB AMR modules. Replicated ranks first
prove exact agreement on root extent, child bounds, refinement ratios, and EB
geometry summaries. The root is split into contiguous y-tiles and each fine
sibling remains one entity. A deterministic greedy scheduler assigns their
64-bit raw, `r`, or `r^2` work to unique owners. The correctness bridge keeps
field replicas but broadcasts each root tile and child state/temperature only
from its owner. Invalid maps or inconsistent work exponents are rejected
collectively before owner-dependent communication, and outputs remain
unchanged on failure.

In `0.111.0`, that bridge executes reactive source integration directly on
owners. The root chemistry field is decomposed over the same y-tiles used by
the distribution, and each fine sibling is advanced as one owner entity.
Covered cells are masked, every owner reactor transaction is accepted
collectively before its state and recovered temperature are broadcast, and a
set-wide fine-to-root average-down restores the composite hierarchy. Inputs
are not published until every entity and the final synchronization succeed,
so a late owner rejection rolls every rank back exactly. The bridge also
rejects rank-inconsistent interval, tolerance, species-width, or
reaction-width controls before owner execution.

In `0.112.0`, the root EB level becomes one exclusive physics entity while
its storage synchronization remains tiled. Weighted StateRedist uses
overlapping, potentially second-order level neighborhoods, so one root owner
preserves the qualified serial redistribution instead of treating tile edges
as artificial boundaries. That owner advances the root once and broadcasts
the state, recovered temperature, and face-centroid fluxes. Each child owner
then performs its own ratio subcycles, accumulates an owner-local flux
register, refluxes the current root candidate, and publishes the corrected
root and child. The root owner performs the final set-wide average-down. Every
stage is collective and the caller's fields are published only after the
complete transaction succeeds.

In `0.113.0`, the same ownership boundary encloses molecular transport. Each
SSPRK2 Euler stage advances the complete root once on its physics owner, then
advances every child on its owner with ratio subcycling, time-interpolated root
exterior data, and a private coarse/fine diffusive flux register. Refluxes are
published in deterministic child order. The root owner applies one set-wide
average-down and, when an interface crosses the embedded boundary, the
qualified global composite conservation closure. The two synchronized Euler
transactions are blended and EOS-recovered on the same owners. Transport
records, boundary data, switches, timestep, and StateRedist controls must agree
bitwise across ranks before any physics starts.

In `0.114.0`, one outer MPI transaction composes owner chemistry, molecular
transport, hydrodynamics, molecular transport, and chemistry over a single
coarse interval. Nested operators work only on private candidate fields. The
caller observes the final hierarchy, transport limiter minimum, and all three
operator counters only after the second chemistry half-step succeeds. A
rejected hydro control therefore rolls back already accepted reaction and
transport prefixes without exposing their local owner counts.

In `0.115.0`, `mpi_amr_eb_sparse_patch_set_2d` separates replicated topology
from numerical payload ownership. Every distribution entry exists as compact
metadata, but a root-tile or child state/temperature allocation exists only on
its owner. Scatter copies owner-authoritative payloads without retaining stale
nonowner fields. Materialization broadcasts each owned tile or child into a
temporary complete patch set and commits it only after collective validation.
This establishes a sparse persistent-storage boundary while legacy physics
still consumes complete temporary arrays.

In `0.116.0`, chemistry consumes the sparse persistent representation
directly. Every root tile and child is reacted only where its owner allocation
exists, using the corresponding EB active mask. Entity acceptance remains
collective and the input sparse object is retained until all local reactors
succeed. Post-reaction fine-to-root average-down currently materializes one
temporary hierarchy on every rank; the result is immediately scattered back
to owner-only payloads before commit.

In `0.117.0`, sparse reactive average-down removes that complete temporary
hierarchy. Each child owner computes the volume-fraction-weighted conserved
state over its coarse footprint and broadcasts only that restriction buffer.
Intersecting root tile owners retain covered cells, recover active-cell
temperature from their local coarse guess, and publish only after every child
is accepted collectively. The child payloads and unrefined root cells remain
owner-local and unchanged.

In `0.118.0`, a sparse outer transaction composes the complete reactive split.
The first chemistry half-step executes and synchronizes directly on sparse
owners. One complete temporary hierarchy then spans both transport half-steps
and the intervening hydro step. The result is scattered back to owners before
the final direct sparse chemistry half-step. Caller state, operator counts,
and the limiter minimum are published only after the entire `R-T-H-T-R`
sequence succeeds.

In `0.119.0`, hydro consumes sparse persistent state without materializing fine
children. Root tiles assemble a level-wide start field because EB StateRedist
still owns overlapping root neighborhoods. The root update and flux arrays are
synchronized, while every child subcycles, accumulates its flux register, and
refluxes only on its owner. Corrected root rows return to sparse tiles and the
existing direct sparse average-down completes the transaction.

In `0.120.0`, SSPRK2 molecular transport also consumes sparse persistent state
without materializing fine children. Each Euler stage assembles only the root
level, advances root transport and StateRedist on its physics owner, and keeps
each child's exterior construction, ratio subcycles, diffusive flux register,
reflux, and state update on that child's owner. A communicator reduction forms
the composite integral needed by cut-interface conservation closure; the
closure corrects only uncovered, unrefined root cells. The two Euler candidates
remain private until the SSPRK2 blend, direct sparse average-down, and collective
validation all succeed.

In `0.121.0`, the sparse outer `R-T-H-T-R` transaction removes its central
replicated compatibility window. It composes the direct sparse chemistry,
SSPRK2 transport, hydro, transport, and chemistry entrypoints on one private
sparse candidate. Each component may use its root-level temporary, but no fine
child payload crosses into a complete replicated hierarchy. The caller state,
all operator counts, and the transport limiter minimum remain unchanged until
every component transaction succeeds.

In `0.122.0`, direct sparse average-down replaces its per-child communicator
broadcast with targeted point-to-point restriction transfer. Each child owner
computes one coarse-footprint buffer, derives the distinct root tile owners
whose y ranges intersect that footprint, and sends only to remote members of
that set. A recipient applies the same conserved restriction and owner-local
EOS recovery as before; unrelated ranks neither allocate nor receive the
buffer. Transfer counts and sparse state publish only after collective
acceptance.

In `0.123.0`, direct sparse hydro removes all-rank root numerical broadcasts.
Sparse root tiles send one packed payload to the root physics owner, which
advances the level-wide root algorithm. That owner sends one packed start,
updated-state, and flux bundle to each distinct remote child owner. Before and
after each remote child reflux, the current root correction moves once in each
direction. After all children succeed, the root owner sends each remote root
tile only its final row band. Full root arrays therefore exist only on the root
owner and ranks that actually own fine children.

In `0.124.0`, direct sparse SSPRK2 transport uses the same targeted root
ownership boundary in both Euler stages. Root tiles gather only to the root
physics owner; each distinct remote child owner receives one packed start,
updated-state, temperature, and diffusive-flux bundle. Cumulative reflux
corrections make one round trip per remote child, and only final row bands
return to remote root tile owners. The final SSPRK2 blend gathers its two root
candidates only to the physics owner and scatters the blended rows. EB-cut
conservation broadcasts only its small conserved boundary-change vector and
applies the uniform correction directly on each locally owned root tile.

In `0.125.0`, timestep selection first consumed the sparse hierarchy through a
targeted root gather. In `0.138.0`, each root tile owner extracts its exact EB
geometry row band and evaluates both hyperbolic and parabolic limits directly
on its exclusively owned state. Each fine child evaluates the same limits only
on its owner and multiplies its local stable step by the refinement ratio to
express a coarse-interval bound. One communicator minimum selects the global
step with no root numerical-field traffic. Fully covered tile bands are
skipped, while an entirely inactive hierarchy rejects rather than publishing
an unbounded interval. Control consensus and field validation happen before
publication, so rejection returns zero dt and zero transfer count without
changing sparse state.

In `0.139.0`, the final sparse SSPRK2 root blend no longer materializes either
candidate on the root physics owner. After the second Euler stage returns its
owned row bands, every root tile owner averages its local interval-start and
second-Euler states, extracts the matching EB geometry band, and recovers
temperature through the mixture EOS. One collective acceptance precedes child
blending and average-down. This removes two root gathers and one row-band
scatter per transport call. The two Euler stages retain their root
gather/scatter boundary because diffusive flux construction and StateRedist
still require neighboring rows.

In `0.140.0`, each sparse transport Euler stage replaces its unconditional
selected-root gather and advance with owner-tiled target-band work. A
target tile owner assembles a six-row guard from its local tiles and direct
row fragments sent by intersecting source owners, then executes the established
EB transport-flux and second-order StateRedist kernels on that band. It routes
only its owned start state, updated state, temperature, x-flux rows, and uniquely
owned y-faces to the root physics owner. That owner still assembles one complete
temporary bundle after tile computation because fine-child exterior
interpolation, flux-register accumulation, and deterministic reflux consume
level-wide arrays. Corrected rows scatter to tile owners before the tile-local
SSPRK2 blend. A target touching a periodic y boundary uses the complete root
band to preserve the serial kernel's cyclic wrap; cyclic finite-halo geometry
is a later optimization.

In `0.141.0`, a periodic y-edge target can replace that complete-root
compatibility band with a boundary-anchored cyclic band. The lower and upper
global boundary remain the temporary geometry's two outer ends, while two
increasing source-row fragments supply the target, its six-row dependency
guard, and one additional seam-isolation row on each side. This prevents the
deliberate internal gap from contaminating a cell on the required stencil.
Absolute EB boundary-centroid y coordinates are shifted into compact-band
coordinates; all other copied EB metrics retain their established local
meaning. If the protected footprint spans the root, the complete-root fallback
remains. The post-compute root bundle for child exterior data and reflux is
unchanged.

In `0.142.0`, an EB flux register allocates correction storage only over the
fine patch expanded by one coarse cell, which contains every coarse/fine flux
mismatch location. Reflux scans only that compact, globally indexed support.
Sparse transport correction traffic expands the patch by one further cell to
cover every cardinal or diagonal recipient of a cut-cell mismatch. A remote
child initializes its full compatibility workspace from the uncorrected root
candidate, receives the latest cumulative values only in that protected
rectangle, performs the established reflux transaction, and returns only the
same rectangle. The root owner merges rectangles in deterministic child order.
The complete start/end/flux input bundle sent once to each distinct child owner
is unchanged and remains the next decomposition boundary.

In `0.126.0`, a public sparse time loop composes that selector with the direct
owner-only `R-T-H-T-R` transaction. The stable interval is recomputed after
every accepted state, the final interval is clipped to the target time, and
time, total step count, minimum accepted dt, operator counts, limiter minimum,
and timestep root-traffic counts publish only after a whole split step commits.
If a later step fails or reaches the configured total-step limit, earlier
committed states and their exact accounting remain visible.

In `0.127.0`, an explicit sparse topology transaction materializes the current
owner fields, applies the established serial multipatch EB average-down,
prolongation, and overlap-retention regrid, recomputes deterministic ownership
for the new ordered children, and scatters the result back to one-copy sparse
storage. Distribution, sparse payloads, and geometry template commit together;
invalid controls or any intermediate failure leave all three unchanged. This
correctness-first compatibility window is confined to regrid events.

In `0.128.0`, temperature-gradient regrid planning is connected to that
transaction and to the public sparse clock. Root tiles gather only to the root
physics owner, which builds the ordered collection and broadcasts compact plan
metadata rather than numerical fields. A caller-supplied geometry builder
reconstructs each planned EB child on all ranks before the existing
serial-compatible regrid commit. When a cadence-triggered regrid is due, the
physics step, distribution, sparse payloads, geometry template, regrid counts,
and root-traffic diagnostics publish atomically. A failed geometry build or
regrid discards the otherwise valid physics candidate and leaves the caller at
the previous accepted time.

In `0.129.0`, the topology transaction removes its replicated numerical-field
window. Old fine owners average down directly to intersecting root-tile owners.
Each distinct new child owner assembles the averaged root from targeted tile
payloads and performs PCM prolongation locally. Same-ratio old/new overlap
rectangles then copy locally or move once from the old child owner to the new
owner before active-cell temperature recovery. Replicated geometry and compact
topology descriptors still define the deterministic transfer schedule. The
candidate distribution, root tiles, children, template, and three transfer
counts publish only after collective validation succeeds.

In `0.130.0`, checkpoint and output consumers gain a root-only sparse
materialization boundary. Each remote root tile and fine child sends one packed
state/temperature payload directly to a caller-selected root. Only that rank
allocates the complete root fields and full patch-set fields; every non-root
output remains unallocated. Root selection and sparse input validity are
collective preconditions, and outputs plus exact sender counts publish only
after the reconstructed hierarchy passes collective validation. The older
all-rank materialization routine remains available for legacy replicated
operators, but checkpoint/output adapters can now avoid rank-multiplied fields.

In `0.131.0`, a dedicated sparse MPI I/O layer owns the writer lifecycle. It
invokes root-only materialization, calls the established formatted multipatch
checkpoint or EB CSV writers only on the selected root, and broadcasts the
writer result before returning. Successful calls publish the sender-local
materialization count; gather rejection, invalid output controls, or root I/O
failure publish zero. Checkpoints retain the serial schema and can be read by
the established reader, while CSV publication writes one root file and one
deterministically named file per child without complete non-root fields.

In `0.132.0`, the inverse boundary reads that checkpoint only on a selected
root. After geometry and topology compatibility checks, the root copies its
owned entities and sends one packed state/temperature payload directly to each
remote root-tile or child owner. Non-root read arrays and patch sets must be
empty on entry. Clock metadata and the sparse candidate publish only after
root agreement, successful read, direct scatter, and collective validation;
failure returns an empty sparse set, zero metadata, and zero traffic.

In `0.133.0`, the restart topology is represented separately from its reactive
fields. Each child descriptor contains only its EB geometry and coarse/fine
patch box; conserved state and temperature are absent from the type. The
root-only reader and direct scatter validate distributions and owner-local
payload shapes against that descriptor, while the selected root alone holds
the complete checkpoint fields. The former full patch-set entrypoints remain
as extraction wrappers. Geometry metadata is still replicated intentionally.

In `0.134.0`, a separate EB patch-tree topology removes the depth limit from
geometry metadata. A root geometry owns a runtime-sized relation sequence;
every relation stores its refinement ratio, flattened ordered children,
parent indices, and parent-to-child offsets. Construction rebuilds every AMR
patch against its actual parent geometry and validates separated siblings. A
dynamic whole-tree replacement is staged as a candidate and commits only after
the complete topology validates; identical plans are explicit no-ops.

In `0.135.0`, a numerical hierarchy mirrors every topology level and patch.
Each node stores its reactive conserved state and recovered temperature while
geometry remains in the topology. Initialization walks parent-to-child and
uses the qualified EB PCM prolongation. Synchronization walks the relations in
reverse and applies the qualified reactive EB average-down to every child.

A dynamic rebuild is one transaction. It first evaluates the old composite
integral and restricts a private old-tree copy to the root. The candidate is
then constructed from that collapsed root. At each new level, children are
prolonged from the already updated parent, same-resolution physical overlap is
retained only when cell and surrounding-face EB metrics agree, and active-cell
temperature is recovered through the NASA7 EOS. A final deepest-to-root
synchronization and complete conserved-vector integral comparison precede the
commit. Invalid plans, incompatible overlap geometry, failed EOS recovery, or
failed conservation leave the accepted topology and every field unchanged.

In `0.136.0`, the replicated MPI-owner hydro path no longer selects one rank
to advance the complete root level. Every root y-tile is extended by at most
six rows on each side, its owner extracts an exact EB geometry band, and the
established reactive EB level kernel advances that bounded band. The guard is
larger than the combined reconstruction, face-centroid interpolation, and
second-order StateRedist dependency radius.

Each tile publishes state, temperature, and x-face fluxes for its owned cell
rows. Y-faces use a unique lower-face ownership rule, with the final tile also
owning the upper physical boundary. Zero-filled rank contributions are summed
to reconstruct the replicated root result. Collective acceptance precedes
assembly, and the caller publishes work counters only after later child
subcycling, reflux, and average-down also commit.

In `0.137.0`, the sparse hydro path uses the same bounded-band numerical
partition without replicating its input. For each target tile, only source
tiles intersecting its six-row band send packed state and temperature row
fragments to that target owner. The owner advances the band and routes its
owned input, updated state, x-flux rows, and uniquely owned y-faces to the root
owner. That rank assembles the complete temporary root bundle required by the
existing child exterior, reflux, and final row scatter stages. Unrelated ranks
never allocate a complete root field. Halo, result, child, correction, and
scatter payloads remain point-to-point, and public advance, work, and traffic
counters remain zero unless the complete sparse transaction commits.

Unsplit transverse prediction, fourth-order StateRedist slopes, periodic ghost
neighborhoods, thermal/viscous/catalytic walls, coarse-to-fine spatial slopes,
same-level diffusive exchange for touching siblings, locally resolved
PeleC-style multilevel redistribution, arbitrary-depth physics recursion,
dynamic root/middle
lifecycle ownership, non-outflow refined boundaries, decomposed sparse root
transport and timestep selection, replica-free regrid overlap transfer,
distributed sparse checkpoint/output, and distributed EB flux registers
remain outside this subsystem.
Dynamic three-level mode changes only the finest patch inside a
fixed middle level and rejects finest removal and siblings.

## Compact sparse child transport context (`0.143.0`)

The sparse transport child phase separates coarse data needed for fine
boundary reconstruction from coarse data needed for reflux. A
`reactive_eb_patch_exterior_context_2d` stores raw start/end conserved state
and temperature samples only on the four fine-patch edges. Physical-boundary
sides retain safe placeholders and are filled from the current fine boundary
state during each substep. Reconstructing from this context uses the same
conserved-state time interpolation and EOS recovery as the complete-root path.

The root physics owner initializes each compact patch-plus-one-cell flux
register and accumulates the coarse interface flux before routing the context
to the child owner. The child performs its ratio subcycles, accumulates fine
fluxes into that register, and returns one fine state/temperature/register
payload. The root owner applies reflux in deterministic child order and sends
back one corrected fine state/temperature payload. A remote child therefore
allocates no complete root start, end, temperature, or x/y-flux array. Hydro
keeps its existing complete child bundle and correction route.

## Compact child-local reactive reflux (`0.144.0`)

The EB reflux kernel accepts a coarse array whose declared lower bounds are
global coarse indices. The array may be the complete root or any rectangle
containing the patch expanded by two coarse cells. Flux mismatch occupies the
patch-plus-one ring; the second cell contains every cardinal or diagonal
recipient of cut-cell redistribution. Reactive temperature recovery scans only
the supplied coarse support and the complete fine patch. Existing complete-
root callers are thin wrappers over the same kernel.

Sparse transport includes the current patch-plus-two coarse state and
temperature in the child context payload. The fine owner subcycles, accumulates
fine interface flux, refluxes locally, retains the corrected fine field, and
returns only corrected coarse support. The root owner merges that support in
ordered child sequence. A remote child transaction therefore contains two
messages per Euler stage and no fine-state round trip.

## Compact coarse interface-flux accumulation (`0.145.0`)

Coarse EB flux-register accumulation accepts x-face and y-face arrays with
explicit global lower bounds. Each rectangle need contain only the active
coarse/fine interface faces: the two vertical face ranges and the two
horizontal face ranges. Validation rejects a missing interface, an out-of-root
face bound, or nonfinite data before modifying the register.

The complete-root entrypoint delegates to this support kernel. Sparse
transport passes patch-local sections of the temporary root flux bundle, so
the consumer no longer requires level-wide array shapes. The temporary bundle
is still assembled on the root physics owner; distributing those interface
fragments directly from root-tile owners is the next ownership boundary.

## Direct root-tile coarse-flux routing (`0.146.0`)

Every sparse root transport tile retains the x-flux rows that correspond to
its owned cell rows. Y-faces follow the established unique ownership rule: a
tile owns its lower face through the face below its upper cell, and the final
tile also owns the upper physical boundary. These retained arrays use global
y lower bounds.

For each child, every intersecting tile owner sends one packed x/y fragment
directly to the child owner. The receiver assembles globally indexed compact
face rectangles and verifies complete coverage before initializing the coarse
flux register. The root physics owner sends only exterior start/end samples
and patch-plus-two state/temperature support; the register is absent from that
message. Child reflux and ordered corrected-support return remain unchanged.

The complete root result and flux bundle still exist on the root physics owner
for exterior extraction, boundary closure, cumulative support merge, and row
scatter. The direct flux route removes one dependency from that boundary but
does not yet distribute state context or corrected support.

## Compact exterior state-context support (`0.147.0`)

Reactive child exterior-context extraction accepts coarse start/end conserved
state and temperature arrays with explicit global lower bounds. A rectangle
containing the patch expanded by one coarse cell covers every nonphysical
coarse/fine boundary sample. Physical sides retain the established placeholder
values and are filled from fine boundary state during subcycling.

The complete-root entrypoint validates its compatibility arrays and delegates
to the support kernel. The unit gate reconstructs exterior arrays from complete
and strictly smaller support contexts and requires bitwise parity. Sparse MPI
continues to call the complete-root wrapper until root-tile state fragments and
cumulative corrected support are routed directly.

## Direct root-tile state/support routing (`0.148.0`)

Each root transport tile retains stage-start state, uncorrected stage-end
state, and current corrected state with their temperatures. A child owner
assembles globally indexed patch-plus-two rectangles directly from every
intersecting tile owner. The start and uncorrected-end views produce the same
four-edge interpolation context as the complete-root path, while the corrected
view carries cumulative reflux changes from earlier children.

Child-local reflux returns one corrected fragment directly to each intersecting
tile owner. Children remain ordered, so overlapping support is visible before
the next child begins. Final corrected root tiles commit locally and the former
root-owner row scatter is absent. Complete temporary root state and flux arrays
remain on the root physics owner for compatibility checks and cut-boundary flux
closure; eliminating that post-compute assembly is a later boundary.

## Owner-local root transport result (`0.149.0`)

Sparse transport Euler stages retain their result only in root-tile state and
flux records. Remote tile owners no longer pack stage-start state, stage-end
state, temperature, or complete owned flux rows for the root physics owner.
The root physics owner therefore allocates no complete transport result or flux
array. Halo input, direct child state/flux support, ordered child correction,
and tile-local final publication keep their established ownership.

Cut-interface conservation needs only the physical root-boundary flux change.
Every tile owner accumulates its left and right x-face contribution. The first
and final y tiles additionally accumulate the lower and upper physical faces.
A communicator sum combines that `nvar` vector on all ranks before the existing
tile-local conservation closure. Hydro and explicit materialization/output
boundaries retain their existing complete-root behavior.

## Compact sparse hydro child context (`0.150.0`)

The root physics owner still assembles the owner-tiled hydro result. It no
longer sends that complete result to each distinct child owner. For every
child, it extracts the established four-edge stage-start/stage-end context,
the current patch-plus-two corrected coarse state and temperature, and only
the x/y face rectangle intersecting the child's coarse boundary. One packed
message transfers those arrays when the child is remote.

The child owner reconstructs the established time-interpolated exterior from
the compact context, accumulates coarse interface flux through the globally
indexed support API, performs ratio subcycling, and refluxes the compact state
support locally. Only corrected patch-plus-two state and temperature return to
the root owner. Child order remains deterministic, so later contexts observe
earlier reflux corrections. Final corrected-root scatter and hydro tile-result
assembly remain subsequent ownership boundaries.

## Direct hydro coarse-flux routing (`0.151.0`)

Each sparse hydro root tile retains the x-flux rows for its cell rows and the
same uniquely owned y-faces used by transport. A child owner assembles its
globally indexed interface rectangles from only the intersecting tile owners.
Local fragments copy directly; remote fragments use one packed x/y message per
tile/child intersection. Complete finite coverage is required before the
coarse register is accumulated.

The tile-to-root hydro result now contains only stage-start and stage-end state
and temperature. The root-to-child hydro context likewise contains only the
four-edge start/end context and current patch-plus-two correction support. The
root physics owner therefore allocates no complete hydro x/y flux array.
Complete root state assembly, ordered corrected-support merge, and final row
scatter remain because the root owner still extracts later child contexts and
publishes the corrected root state.

## Owner-local hydro result and direct state/support routing (`0.152.0`)

Hydro root tiles retain stage-start, uncorrected stage-end, and current
corrected state and temperature beside their owner-local flux records. For each
child, intersecting tile owners send only patch-plus-two row fragments to the
child owner. That owner assembles globally indexed start/end/corrected support,
checks complete coverage, and extracts the established four-edge interpolation
context locally.

After child-local subcycling and reflux, corrected support returns directly to
each intersecting root tile owner before the next child is assembled. Final
corrected root rows publish locally, so sparse hydro has no complete root state,
temperature, or flux result, no remote tile-result message, and no final root
scatter. Finite-band halos, deterministic child order, direct flux fragments,
and hierarchy-wide average-down remain unchanged. Superseded private root-
bundle, root-context, correction, and scatter communication helpers are absent.

## Arbitrary-depth reactive EB patch-tree timestep (`0.153.0`)

The single-node active-cell CFL calculation now lives below both the runnable
driver and AMR hierarchy layers. The existing driver entrypoint remains a thin
compatibility wrapper, while the reactive EB patch tree calls the same kernel
for every root and child node without materializing another hierarchy.

Tree traversal carries the cumulative product of relation refinement ratios.
Each node-local interval is multiplied by that product before entering the
root-time minimum, matching the number of temporal subcycles from that node to
the root. Fully covered nodes impose no stability constraint and are skipped;
an entirely inactive tree rejects. Invalid trees, species layouts, CFL
controls, active-node states, or an overflowing cumulative scale reject with
zero timestep. The accepted hierarchy is read-only throughout selection.
Hydro, chemistry, transport, public clock ownership, and MPI distribution
remain separate arbitrary-depth operations.

## Arbitrary-depth reactive EB patch-tree hydro (`0.154.0`)

Hydrodynamics is a transaction over a private numerical-tree candidate. One
recursive invocation advances one node for its supplied interval and retains
the node's start state, uncorrected end state, and x/y EB fluxes. If the node
has children, each child receives time-interpolated exterior state from those
two parent endpoints and advances recursively for exactly the relation ratio
substeps.

Every parent/child pair owns an independent EB flux register. Coarse flux is
accumulated once for the parent interval; child flux is accumulated after every
recursive substep. Children reflux in topology order, then average down into
their actual parent. Each refined subtree compares its before/after composite
integral against flux through the parent node's outer boundary. Any remaining
density, total-energy, and species residual is distributed over active,
unrefined parent cells with EOS recovery and mass/species closure validation.

After the root recursion succeeds, one deepest-first synchronization restores
all coarse representations before the candidate commits. A failed level
advance, exterior fill, register operation, reflux, conservation closure, EOS
recovery, or final validation leaves the accepted tree unchanged and returns
zero per-level advance counts. Chemistry, molecular transport, a public clock,
dynamic tags, checkpoint I/O, and MPI ownership remain separate.

## Arbitrary-depth reactive EB patch-tree chemistry (`0.155.0`)

The numerical tree now owns a shared active-cell chemistry traversal. Every
runtime patch obtains its own geometry from the topology, masks covered EB
cells, and calls the established 2D constant-volume chemistry integrator once
per requested reaction interval. The standalone operation synchronizes a
private candidate deepest first and commits per-level patch-call counts only
with the complete tree.

The tree Strang entrypoint applies chemistry for `dt/2` on every node, invokes
the recursive hydro transaction for `dt`, applies chemistry for `dt/2` again,
and performs final deepest-first synchronization. Both chemistry and hydro
operate on the same private candidate. A rejection in the first reaction
stage, recursive hydro, second reaction stage, EOS recovery, synchronization,
or validation publishes neither fields nor counters. Molecular transport, a
public clock, dynamic tags, checkpoint I/O, and MPI ownership remain separate.

## Arbitrary-depth reactive EB patch-tree transport (`0.156.0`)

One recursive transport Euler call advances one runtime node, retaining the
node start and Euler-end fields for child-time interpolation. Every child takes
the relation refinement-ratio subcycles, owns an independent diffusive EB flux
register, refluxes in deterministic topology order, and averages down into its
actual parent. Each refined subtree closes density, total energy, and species
against the parent's outer diffusive flux before the final deepest-first
synchronization.

The public SSPRK2 operation runs that complete recursive Euler transaction
twice on a private tree, blends every node with its accepted stage-zero state,
recovers active temperatures through the EOS, and synchronizes again. State,
temperature, the minimum positivity-limiter theta, and optional per-level node
counts publish only after the final candidate validates. A combined
`R-T-H-T-R` transaction, public clock, dynamic tags, checkpoint I/O, and MPI
ownership remain separate.

## Arbitrary-depth reactive EB patch-tree full physics (`0.157.0`)

The full-physics entrypoint owns one private numerical-tree candidate across
active-cell chemistry, recursive SSPRK2 transport, and recursive
hydrodynamics. It applies `R(dt/2)`, `T(dt/2)`, `H(dt)`, `T(dt/2)`, and
`R(dt/2)` in that order, reusing the qualified standalone tree operations.

Chemistry patch calls, transport Euler-node calls, hydro node calls, and the
minimum transport limiter theta accumulate privately. Final deepest-first
synchronization and complete tree validation precede one atomic publication.
Any rejection after an earlier valid physics prefix therefore preserves the
accepted hierarchy and returns zero public counts plus theta one. Dynamic
tags, checkpoint I/O, and MPI ownership remain separate.

## Arbitrary-depth reactive EB patch-tree time loop (`0.158.0`)

The public clock recomputes both active-cell hyperbolic and explicit mixture
transport limits on every runtime node. Each local interval is multiplied by
the cumulative ancestor refinement product before the global tree minimum is
clipped to the remaining target time.

Every accepted interval runs the full `R-T-H-T-R` operation on a private tree
candidate. The candidate tree, time, total step count, minimum accepted
interval, limiter minimum, and accumulated per-level physics counts publish
together. A rejected first step changes nothing; reaching the caller's step
limit after prior success retains exactly that committed prefix. Dynamic tags,
checkpoint I/O, and MPI ownership remain separate.

## MPI arbitrary-depth EB patch-tree ownership (`0.159.0`)

The first distributed layer retains the replicated numerical tree but assigns
each runtime node to one deterministic owner. Greedy placement follows
topology order and minimizes accumulated rank work; a configurable exponent
weights deeper nodes by their cumulative subcycle product. Per-rank cell,
entity, and weighted-work totals are derived from the owner map and validated
against it.

Before publication, every rank must agree on the topology geometry, relation
ratios, and weighting control. Each owner broadcasts its candidate state and
temperature into a private replicated tree, and all ranks commit only after
the complete candidate validates. Rank-local invalid input or inconsistent
controls reject collectively with zero publication accounting. Sparse field
storage, direct owner migration, and owner-local physics remain separate.

## MPI sparse arbitrary-depth EB patch-tree storage (`0.160.0`)

The sparse numerical tree retains the complete topology on every rank but
allocates conserved state and temperature only for locally owned nodes.
Initialization copies each node from the accepted replicated tree to its
owner; nonowners retain no field allocation. A deliberate materialization
boundary broadcasts owner-authoritative fields into one private complete tree
and publishes it only after collective validation.

Ownership changes allocate a second sparse layout from the new map. Nodes
whose owner is unchanged copy locally, while changed nodes send state and
temperature directly from the old owner to the new owner. The old sparse tree
is replaced only after every rank validates the complete candidate. Invalid or
inconsistent owner metadata rejects before transfer and preserves the accepted
sparse tree exactly. Distributed timestep reduction and owner-local recursive
physics remain separate.

## MPI owner-local arbitrary-depth EB patch-tree timestep (`0.161.0`)

The combined hyperbolic and explicit-transport stability scan now consumes the
sparse tree directly. Each rank visits only its owned active nodes, evaluates
the qualified EB hydro and transport limits, and multiplies each node interval
by its cumulative ancestor refinement product. Ranks with no active ownership
contribute a neutral huge value.

The communicator minimum is published only after collective distribution,
field, species-layout, CFL, and transport-control checks. The sum of evaluated
nodes must be nonzero, and the final root interval must be finite, positive,
and non-huge. No replicated numerical tree is constructed. Recursive hydro,
transport, and chemistry execution remain separate owner-local boundaries.

## MPI owner-local arbitrary-depth EB patch-tree chemistry (`0.162.0`)

Chemistry now advances one private sparse candidate. Every rank traverses the
same node order, but only the assigned owner runs the active-cell reactor and
temperature recovery. A collective accept boundary follows every node so a
rank-local failure cannot expose an accepted prefix or desynchronize later
communication.

Hierarchy synchronization then traverses relations deepest-first. When child
and parent share an owner, average-down is local. Otherwise the child owner
sends conserved state directly to the parent owner, which applies the same EB
restriction kernel as the serial tree. Every child has a collective accept
boundary, and the sparse candidate commits only after final all-rank
validation. Recursive hydro and transport remain separate owner-local work.

## MPI sparse arbitrary-depth EB composite integrals (`0.163.0`)

The full-tree integral is a root-subtree wrapper. The subtree operation first
establishes communicator agreement on the replicated ownership/topology,
sparse layout, selected level and patch, and conserved-component extent. A
recursive local walk builds each node's direct-child refined mask, integrates
only unrefined cells when that node is owned locally, and always descends
through the replicated relation graph.

One `MPI_SUM` combines the rank-local conserved vectors and a second sum checks
the number of contributing owner nodes. No numerical node is allocated or
materialized on a nonowner. Invalid selectors, rank-dependent selectors,
nonfinite results, or an empty contributing set publish a zero integral and
zero optional local-node count. The same subtree boundary can therefore be
used before and after owner-local reflux and cut-interface closure.

## MPI owner-local arbitrary-depth EB patch-tree hydro (`0.164.0`)

Every rank follows the serial depth-first subcycle schedule while only the
selected node owner executes its EB level update. Before an owner boundary is
crossed, the parent owner extracts the compact start/end four-edge exterior
context and sends it once to the child owner. Each fine substep returns its
x/y flux vector directly to the parent owner, which retains and consumes that
edge's coarse/fine flux register.

Reflux keeps the parent field on its owner: a distinct child owner sends its
current node, the parent applies the serial reactive reflux kernel, and the
corrected child returns. Average-down then sends the corrected child once more
to the parent in deterministic child order. Shared-owner edges execute every
operation locally. Each refined node measures its subtree integral before and
after the operation and applies the established unrefined-parent conservation
closure on the parent owner. These internal reductions reuse the already
validated topology and avoid repeating public metadata consensus at every
subcycle.

All fields, per-level advance counts, and grouped direct-transfer counts remain
inside one private sparse candidate until the complete recursive root operation
validates collectively. A control mismatch or later owner failure therefore
publishes zero accounting and preserves the accepted sparse tree exactly.

## MPI owner-local arbitrary-depth EB patch-tree transport (`0.165.0`)

Each SSPRK2 Euler stage follows the recursive hydro ownership schedule, but the
node owner evaluates molecular-transport fluxes, the conservative RHS, and
StateRedist. Parent start/end exterior context crosses a distinct-owner edge
once per node invocation, and each fine substep returns its diffusive fluxes
directly to the parent-owner register. Reflux and ordered average-down reuse the
same direct child-state routes and subtree conservation closure as hydro.

The second Euler stage advances a private copy of the first. Every owner then
blends its own start and second-stage fields and recovers temperature locally.
One final deepest-first restriction makes covered parent cells authoritative;
no complete numerical node or tree is constructed on a nonowner. The public
minimum limiter is an `MPI_MIN` over owner-local values.

Boundary data, transport flags, interval, redistribution controls, species
layout, topology, ownership, and sparse fields must agree before advancement.
The two Euler stages, final blend, hierarchy synchronization, limiter minimum,
per-level advances, and grouped direct-transfer counts publish only after the
complete sparse candidate validates. Full-physics composition and the public
sparse clock remain separate transactions.

## MPI owner-local arbitrary-depth EB full physics (`0.166.0`)

The sparse `R-T-H-T-R` entrypoint owns one private numerical-tree candidate.
It applies an optional chemistry half-step, an SSPRK2 transport half-step, the
complete recursive hydro interval, a second transport half-step, and a final
optional chemistry half-step. Every stage calls the qualified owner-local
operator directly; no complete tree is materialized between stages.

An outer communicator preflight makes the timestep, tolerances, redistribution
controls, physics flags, and species/mechanism/transport extents identical
before any rank branches on optional physics. Inner stage consensus continues
to validate boundary data and hydro scheme strings. A rejected later stage
therefore discards earlier valid prefixes with the accepted sparse tree and
all public diagnostics unchanged.

Per-level chemistry, transport-Euler, and hydro advances accumulate separately.
Restriction, transport-route, and hydro-route transfers likewise remain
separate so the topology/owner map predicts each category exactly. Both
transport limiter minima reduce to one public value. Fields and every counter
commit only after final sparse validation. Target-time clock ownership remains
separate.

## MPI owner-local arbitrary-depth EB target-time clock (`0.167.0`)

The public sparse clock first establishes exact communicator agreement on the
accepted time and step, target time, step ceiling, CFL values, solver controls,
physics flags, and data extents. It then evaluates the qualified owner-local
hydro/transport timestep before every attempted step and clips the result to
the remaining target interval.

Each step advances a private sparse candidate through the owner-local
`R-T-H-T-R` transaction. Only after that candidate and all category counters
validate does the clock commit fields, time, step count, minimum accepted dt,
minimum transport limiter, timestep-node evaluations, per-level advances, and
operator-specific transfer counts. Successful completion assigns the requested
target time exactly.

If the step ceiling is reached or a later timestep/physics operation fails,
the already committed prefix remains authoritative with matching diagnostics;
the uncommitted step is discarded. An initial rank-dependent clock control
rejects before timestep evaluation with neutral outputs. Dynamic tagging,
checkpoint/restart, and output for this arbitrary-depth tree remain separate.

## Serial arbitrary-depth EB temperature-tagged rebuild (`0.168.0`)

The serial planner first copies and deepest-first synchronizes the accepted EB
tree. It then visits every parent at one prospective relation, applies the
existing normalized temperature-gradient tagger and disconnected-component
clusterer, and records children in parent-major deterministic order. Each
rectangle is converted to EB geometry by a caller-supplied builder so the
planner remains independent of a particular level-set representation.

After one relation is accepted, a temporary topology and PCM-prolongated field
tree provide the prospective parent temperatures for the next relation. This
continues until the level ceiling, no tags, or no taggable parent remains. The
accepted solution is never modified during planning.

The public regrid wrapper passes the complete plan to the established
overlap-preserving transactional rebuild. Identical topology is a no-op;
empty plans collapse to the root; changed plans retain geometrically matching
same-resolution overlap and conservatively initialize the remainder. Any tag,
geometry, EOS, topology, or conservation failure leaves the accepted tree
unchanged. MPI owner-local planning and topology-changing migration remain a
separate transaction.

## MPI owner-local arbitrary-depth EB temperature-tagged rebuild (`0.169.0`)

The sparse planner copies topology and geometry but leaves every numerical
field on its current owner. For each prospective parent, only that owner runs
the temperature tagger. Integer reductions expose compact tag counts and
bounds so every rank reconstructs the same deterministic parent-major plan and
invokes the caller's geometry builder in the same order.

The candidate topology receives a new deterministic work-weighted owner map.
PCM initialization sends parent state directly to each new child owner.
Geometrically identical old/new rectangles retain same-resolution cells by
direct old-owner to new-owner transfers; no complete field tree is gathered.
A deepest-first direct restriction then closes every parent/child relation.

Topology, distribution, state, temperature, tagged/transferred-cell counts,
and restriction traffic commit together only after topology/geometry checks,
EOS recovery, sparse validation, and a composite-integral test succeed on all
ranks. An unchanged topology is a field-exact no-op. Empty tag plans collapse
the tree to its synchronized root. Arbitrary-depth checkpoint/restart and
composite output remain separate lifecycle work.

## Serial arbitrary-depth EB patch-tree checkpoint (`0.170.0`)

The checkpoint is a distinct versioned formatted stream. Its header fixes the
species order, conserved-state extent, and level count. It then stores the root
EB geometry and every ordered relation: refinement ratio, parent index,
coarse-cell rectangle, and complete child EB geometry. Geometry records include
cell volumes and centroids, face apertures and centroids, embedded-boundary
lengths, centroids, normals, normal integrals, and cell classifications.

Lifecycle time, minimum accepted timestep, step count, and regrid count precede
the level-major, patch-major numerical fields. A terminal marker detects
truncation. Reading constructs topology and fields in a private candidate,
recovers temperature from conserved state with the selected species database,
and publishes only after complete structural and thermodynamic validation.

Level, patch, and geometry-cell limits are checked before allocation. Schema,
species-order, depth, topology, dimension, finite-value, EOS, or end-marker
failure returns an empty tree and zero metadata. Sparse MPI checkpoint I/O and
rank-neutral restart remain separate lifecycle work.

## Sparse MPI arbitrary-depth EB checkpoint/restart (`0.171.0`)

The write boundary validates communicator-wide root, clock metadata, counters,
and species order. Every node already owned by the selected I/O root is copied
locally; each remaining owner sends its state and temperature once directly to
that root. Only the root constructs the complete serial tree and writes the
qualified `0.170.0` format.

On restart, only the selected root reads and validates the file. It broadcasts
level/relation rectangles and complete EB geometry, not numerical fields. Each
rank reconstructs the same topology and computes a fresh work-weighted owner
map for the current communicator and requested subcycle exponent. The root then
sends each non-root-owned node directly to its new owner.

One sender-side entity transfer is counted per root/owner difference in either
direction. A rank-dependent root, depth, exponent, lifecycle value, or species
order rejects before file I/O or field traffic. Read, topology broadcast,
distribution construction, sparse scatter, and final validation must all
succeed before any public output becomes nonneutral. Composite output remains
separate lifecycle work.

## Arbitrary-depth EB composite output (`0.172.0`)

The serial writer traverses every level and patch in deterministic order and
constructs a coarse-cell mask from that patch's direct children. It writes only
unmasked cells, so a parent cell replaced by any finer child is omitted while
every finest available cell is emitted once. Each row identifies its level,
patch, local indices, spacing, physical center, EB volume fraction and boundary
metrics, conserved density and total energy, recovered primitive state,
temperature, and ordered species mass fractions.

The sparse MPI adapter validates collective agreement on writer root, time,
and species order, then uses the existing direct node gather. Only the selected
root materializes the complete numerical tree and opens the CSV. Every remote
node contributes one entity transfer; non-root ranks retain only their owned
fields. The root broadcasts the final write status, and transfer counts remain
neutral when control, gather, thermodynamic conversion, or file output fails.

## Runnable serial arbitrary-depth EB application (`0.173.0`)

`pelef_reactive_eb_patch_tree_2d` is a separate public executable so the
legacy single-patch, sibling-multipatch, and fixed three-level application
contracts remain unchanged. It reuses their reactive-flow, embedded-boundary,
and AMR namelists and adds only `patch_tree_maximum_levels`.

A fresh run constructs the configured root EB geometry and reactive field,
initializes a root-only tree, and optionally applies recursive temperature-tag
planning before the first step. Each committed root step selects the minimum
hydro/transport limit over all nodes, advances the qualified `R-T-H-T-R`
transaction, then applies scheduled topology rebuild and checkpoint output.
Restart delegates to the self-describing tree reader. Completion, checkpoint
stop, or restart all use the same single composite CSV writer.

The driver owns time, step, regrid, minimum-dt, and transport-limiter state.
Numerical topology and fields remain owned by the patch-tree core; geometry
construction is an internal callback using the configured plane or circle over
each child region. Existing fixed-depth application modes do not call this
path.

## Public patch-tree checkpoint/restart lifecycle (`0.174.0`)

The public serial application now has a process-boundary qualification. The
checkpoint-stop process performs initial recursive tagging, commits one
full-physics root step, applies its scheduled recursive regrid, writes the
self-describing tree, and exits after publishing a composite CSV. A second
process reconstructs the hierarchy and fields from that file and resumes the
same global root-step cadence.

The checkpoint owns numerical state: tree geometry and relations, every node
field, time, committed root-step count, regrid count, and minimum accepted
timestep. Continuation controls such as final time, CFL, physics switches, and
tagging thresholds remain explicit in the restart input. The parity gate uses
identical continuation controls and compares the restarted result with an
uninterrupted reference by stable `(level, patch, i, j)` identity.

## Public sparse-MPI patch-tree application (`0.175.0`)

`pelef_mpi_reactive_eb_patch_tree_2d` composes the qualified sparse ownership,
tagging, timestep, full-physics, regrid, integral, checkpoint/restart, and CSV
APIs behind the established reactive 2D/EB/AMR namelists. A configurable work
exponent weights deeper nodes during deterministic owner assignment.

A fresh run constructs only the replicated root field, converts it to sparse
ownership, releases the replicated state, and performs recursive initial
tagging owner-locally. Thereafter numerical node fields remain allocated only
on their owners. Checkpoint and output gather only to selected root zero;
restart reads there and scatters directly under the current rank count. The
temporary replicated root initialization remains an explicit startup boundary.

## Public sparse-MPI cross-rank restart (`0.176.0`)

The installed sparse-MPI application now has a process- and ownership-boundary
qualification. An uninterrupted one-rank process constructs and advances the
four-level hierarchy with depth-weighted ownership. A separate two-rank
process uses uniform node weighting, writes the self-describing checkpoint
after one committed root step, publishes its intermediate composite, and
stops. Independent four- and eight-rank processes read that checkpoint,
recompute depth-weighted ownership for their communicator, scatter fields
directly to those owners, and resume the global root-step cadence.

No checkpoint owner map is authoritative across the boundary. Geometry,
relations, fields, time, committed-step count, regrid count, and minimum
accepted timestep come from the checkpoint; communicator size and the MPI
work exponent come from the restart process. Final composite comparison by
stable `(level, patch, i, j)` identity therefore covers both rank-count and
ownership-policy redistribution without permitting a replicated numerical
child tree.

## Owner-local public sparse-MPI startup (`0.177.0`)

Fresh application startup now constructs the geometry-only root topology and
its deterministic distribution before allocating numerical fields. Exactly
the owner of root node `(level=0, patch=1)` calls the established reactive 2D
initializer. Non-owners retain unallocated root state and temperature
variables throughout startup.

The root-only sparse initializer collectively verifies a one-level, one-node
topology, communicator-consistent state width, owner-only input allocation,
exact field shapes, finite values, and positive temperature. It then transfers
the owner's allocatable state and temperature directly into the sparse node
with `move_alloc`; no field copy or numerical broadcast occurs. The public
driver requires exactly one initializer rank and requires both source arrays
to be unallocated after the transfer before recursive owner-local tagging can
begin. Replicated EB geometry and tree relations remain intentional compact
metadata needed for deterministic planning and routing.

## Public patch-tree checkpoint fingerprint (`0.178.0`)

Public serial and sparse-MPI checkpoint writes use schema 2 and place a
structured compatibility fingerprint after the ordered species header. It
records root mesh/domain, EB geometry parameters, hierarchy and refinement
controls, numerical method names, chemistry/transport switches and tolerances,
StateRedist controls, and dynamic tagging/regrid controls. Restart compares
integer and character fields exactly and round-trip real fields within a small
machine-precision bound before reading geometry or numerical payloads.

Evolved topology, fields, clock, and lifecycle counters remain checkpoint
state. Final time, maximum steps, output paths, checkpoint cadence, MPI rank
count, and MPI work exponent remain restart-mutable controls. Thus the existing
two-to-four/eight-rank redistribution remains valid, while changing a physics
control such as CFL rejects the file transactionally. Low-level verification
callers without a fingerprint retain the isolated schema-1 compatibility API.

## Interface-local multilevel EB conservation closure (`0.179.0`)

Every direct child rectangle contributes a coarse-side interface support. For
each coarse cell immediately outside that rectangle, the closure marks the
active, unrefined cells in its clipped three-by-three neighborhood. The union
over siblings is the only admissible recipient set for density, total-energy,
and species residuals after reflux and average-down. Physical-boundary sides
without a coarse neighbor contribute no support.

The recipient set is derived from replicated geometry and topology, so serial,
fixed-depth, multipatch, arbitrary-depth, and sparse-MPI paths select the same
cells. Sparse MPI applies the correction only on the owning parent node or root
tile. Recipient fluid volume normalizes the correction; EOS recovery and the
existing final composite-integral check remain transactional. This removes the
former parent-wide perturbation but does not claim bitwise equivalence to
AMReX's per-neighborhood `MLStateRedistribute` transfer bookkeeping.

## Embedded-wall molecular heat and momentum transfer (`0.180.0`)

`reactive_boundary_set_2d` owns one validated embedded-wall record in addition
to its four Cartesian domain faces. Its default is a stationary adiabatic slip
wall with zero species flux, preserving every earlier EB result. A caller may
select an isothermal wall temperature, no-slip velocity, or both through that
record without changing any transport stepping interface.

For each cut cell, `eb_reactive_transport_2d_mod` recovers the cell primitive
state and mixture transport coefficients, measures the centroid-to-wall
distance along the solid-to-fluid normal, and forms one wall-normal flux. The
wall length converts that flux to an extensive contribution and the cut-cell
fluid volume converts it to the local right-hand side. Existing StateRedist,
EOS recovery, AMR reflux, sparse ownership, collective validation, and rollback
remain downstream of the same source.

The MPI control-consensus paths compare the embedded-wall strings, temperature,
velocity, and allocated boundary vectors in addition to the four domain faces.
The low-level boundary-set API is qualified across every current transport
path. Namelist exposure and checkpoint fingerprinting of nondefault wall values
remain a separate public-application lifecycle milestone.

## Public single-level embedded-wall controls (`0.181.0`)

The `&embedded_boundary` namelist now owns the wall kind, thermal mode,
temperature, and three-component velocity used by the public single-level EB
application. Configuration validation couples an isothermal selection to
enabled thermal conduction, no-slip to enabled viscosity, and nonzero wall
velocity to no-slip. The boundary builder applies these values transactionally
only after its domain faces and embedded-wall storage are valid.

The single-level application is checkpoint-free, so this exposes the wall
physics without creating an untracked restart dependency. The AMR public
driver explicitly rejects an active isothermal or no-slip embedded-wall config
at preflight while its checkpoint/fingerprint formats remain unchanged. The
low-level AMR and MPI boundary-set APIs qualified in `0.180.0` remain available
to library callers.

## Restart-safe AMR embedded-wall controls (`0.182.0`)

All public AMR drivers now construct their domain and embedded-wall boundary
records through the same configured builder as the single-level driver. The
two-level, sibling-patch, three-level, arbitrary-depth, and sparse-MPI paths
therefore receive identical wall kind, thermal mode, temperature, and velocity
values before any state allocation or advancement.

The fixed-depth formatted checkpoint schemas advance to version 2 and store
those wall values together with the transport enable flag, individual
viscosity, conduction, diffusion, and barodiffusion flags, and transport CFL.
The serial/sparse arbitrary-depth fingerprint advances to schema 3 and compares
the same controls. A restart mismatch returns transactionally with a neutral
clock and no candidate solution. Earlier schemas are rejected rather than
silently assuming the new defaults.

## EB-safe limited-linear AMR prolongation (`0.183.0`)

`prolong_reactive_eb_patch_linear_2d` computes component-wise monotonized-
central slopes from active regular coarse neighbors. Cartesian fine-child
offsets have zero parent mean, so every accepted regular-parent interpolation
restricts to its source conserved state without a correction pass. Conserved
states, rather than temperatures or primitive variables, are interpolated;
temperature is recovered independently through the configured EOS.

A parent receives linear slopes only when it and all of its fine children are
regular. Cut, covered, or topology-mismatched parents use the established PCM
state. If any linearly reconstructed child is outside the EOS-admissible set,
the complete parent is retried with PCM before publication. Invalid inputs
leave both output arrays neutral. At the `0.183.0` boundary, public regrid
orchestration continued to use PCM pending a separate input lifecycle.

## Fixed-depth public prolongation selection (`0.184.0`)

`reactive_eb_amr_2d_config` now carries `prolongation_method`, read from the
public `&eb_amr` namelist and restricted to `pcm` or `linear`. One dispatcher
owns that method boundary, while the existing low-level PCM and limited-linear
kernels retain their separate numerical contracts. Static initialization and
dynamic replacement of a two-level fine patch, sibling-patch set construction
and replacement, and both three-level coarse-to-middle and middle-to-finest
initializations all pass the selected method explicitly.

The default remains PCM for backward compatibility. The public hot-wall AMR
transport regression selects linear and therefore exercises the installed
application path rather than only a library call. Fixed-depth checkpoint
formats and the arbitrary-depth checkpoint fingerprint do not yet store the
method. Configuration preflight consequently rejects linear whenever a
fixed-depth checkpoint or restart path is active, and arbitrary-depth serial
and sparse-MPI orchestration rejects non-PCM before allocating a candidate.

## Restart-safe arbitrary-depth prolongation selection (`0.185.0`)

The single-patch, sibling-patch, static three-level, and dynamic three-level
formatted checkpoint schemas advance to version 3. Each writes the selected
prolongation method beside the reconstruction controls and compares it before
reading topology or numerical fields. A mismatch leaves the result neutral.

The serial/sparse patch-tree fingerprint advances to schema 4 and owns the
same method string. Serial iterative tag planning, final tree rebuilding, and
sparse owner-local candidate planning and rebuilding now pass the selection to
the shared prolongation dispatcher. Sparse MPI encodes `pcm` and `linear` as a
collective control, rejects invalid or rank-inconsistent selections, and only
the parent owner constructs each new child before existing direct routing.
The established overlap retention, EOS recovery, synchronization,
conservation checks, and atomic publication remain downstream.
