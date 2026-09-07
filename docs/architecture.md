# PeleF architecture

## Executable split

PeleF exposes seventeen default serial verification drivers, up to nine
configure-time selected-mechanism serial drivers, and eleven optional MPI
drivers over shared numerical and physical-property modules.

```text
pelef
  └─ one-dimensional constant-gamma PCM / PLM / PeleC-style paths

pelef2d
  └─ two-dimensional periodic constant-gamma CTU-style Euler path

pelef3d
  └─ three-dimensional periodic constant-gamma PCM / SSPRK2 Euler path

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

pelef_reactive_3d
  └─ periodic NASA7 multispecies PCM hydro, chemistry, and transport

pelef_reactive_eb_3d
  └─ stabilized planar NASA7 EB hydro with FluxRedist or StateRedist

pelef_amr_reactive_3d
  └─ static two-level NASA7 hydro with subcycling, reflux, and average-down

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

## Conservative limited-linear cut-parent prolongation (`0.186.0`)

The shared limited-linear kernel now treats an EB-cut parent as a fluid-volume
control volume instead of forcing it to PCM. It measures directional
differences between active coarse fluid centroids, uses MC limiting when both
sides exist, and retains a one-sided derivative next to covered geometry. A
component-wise envelope limiter bounds reconstructed children by the active
3-by-3 coarse neighborhood.

Fine fluid-centroid offsets generally do not have zero mean within a cut
parent. The kernel removes their fine-volume-fraction-weighted mean before
reconstruction, so EB average-down recovers the source conserved state to
roundoff. Covered parents and regular parents with topology-mismatched child
blocks retain PCM. EOS recovery and the existing parent-local PCM retry remain
the final transactional acceptance boundary. Because fixed-depth, patch-tree,
serial, and sparse-MPI lifecycles already share this dispatcher, they inherit
the cut-parent behavior without a new checkpoint identity field.

The end-to-end reactive qualification also closes the interaction with
second-order weighted StateRedist. If its reconstructed conserved state fails
EOS recovery, the provisional state is redistributed again with order zero.
Both candidates use the same neighborhood partition, so the retry retains the
componentwise volume-conservation contract and the transaction commits only an
EOS-admissible whole-level state.

## Multidimensional cut-parent prolongation (`0.187.0`)

The cut-parent gradient now comes from one least-squares system over active
fluid-volume centroids in the coarse 3-by-3 neighborhood. Axial neighbors must
share an open face. A diagonal neighbor must have at least one two-face path
through an active intermediate cell, which prevents the fit from crossing a
covered corner. A full-rank normal matrix supplies both gradient components;
a rank-one matrix uses the minimum-norm gradient along its only resolved
direction, and an empty stencil remains constant.

The coarse active-neighbor envelope limits predictions at connected coarse
centroids, then the existing fine-child envelope limits the actual
reconstruction. Both limiters preserve an exact affine field when its child
predictions remain inside that envelope, including the qualified
interface-tangential field. The existing fine-volume-fraction-weighted
zero-mean offset, EOS recovery, parent-local PCM retry, and shared lifecycle
dispatcher remain unchanged.

## Rank-recovering cut-parent prolongation stencil (`0.188.0`)

The cut-parent fit first retains the compact connected 3-by-3 stencil. When
that normal matrix is rank deficient, the kernel rebuilds the complete system
over the parent-centered 5-by-5 box. A bounded face-connectivity flood fill
admits only active cells reachable from the parent through open Cartesian
faces inside that box; covered cells and disconnected fluid components cannot
enter the fit.

The selected stencil supplies both the normal equations and the component
envelope. A full-rank grown system restores both affine-gradient components;
only a still-rank-deficient grown system uses the minimum-norm rank-one
fallback. The fine-volume-weighted zero-mean correction, fine-child limiter,
EOS recovery, PCM retry, and shared lifecycle dispatcher remain unchanged.

## Transactional fixed three-level parent regrid (`0.189.0`)

Moving the root-to-middle rectangle changes the coordinate system that owns
the finest patch. The fixed-depth library therefore treats the complete
replacement as one candidate. It restricts the old finest state into the old
middle, regrids the synchronized middle against the root, and retains any
same-resolution middle overlap. It then tags the rebuilt middle only inside
the established two-cell safety margin and prolongs a replacement finest
patch.

The original root, middle, finest, geometries, and patch descriptors remain
untouched until both transfers validate and the before/after three-level
composite conserved integrals agree. A missing valid interior finest plan is a
transaction failure because the fixed three-level representation cannot
publish a temporarily inactive finest level. The public schedule and
fixed-depth checkpoint schema do not yet activate this parent transaction.

## Public fixed three-level parent lifecycle (`0.190.0`)

`dynamic_parent_regridding` opts the fixed three-level application into a
parent-first topology schedule. Initialization and each accepted regrid
interval first call the complete parent transaction. A changed parent already
contains a rebuilt finest patch and increments the counter once. An unchanged
parent falls through to the established finest-only transaction, preserving
the old behavior for a stationary root plan.

Dynamic three-level checkpoint schema 4 records the policy flag and both
actual patch descriptors. Restart validates the stored parent against the root
domain, derives the middle dimensions from that descriptor, validates the
finest two-cell margin in those derived dimensions, and reconstructs both EB
geometries in private candidates. State, temperature, geometry, counters, and
patch descriptors publish together only after EOS recovery and the terminal
marker succeed. Inputs that omit the flag retain a fixed configured parent.

## Arbitrary-depth outflow-boundary children (`0.191.0`)

The serial and sparse-MPI public patch-tree applications now qualify the
domain-inclusive rectangles already produced by the shared temperature-tag
planner. A recursively tagged x-upper hotspot builds four populated levels
whose child rectangles all meet the same physical boundary.

Each physical child side obtains its exterior state by copying the current
fine boundary cell, matching the established outflow contract. Coarse-time
interpolation remains active on the other sides. Flux-register accumulation
and reflux remain absent on the physical side while the other coarse/fine
interfaces retain their normal conservative synchronization. Sparse ownership
uses the same compact context and direct owner routes; no complete field is
introduced on a nonowner rank.

## Boundary-touching patch-tree restart (`0.192.0`)

The arbitrary-depth checkpoint already stores domain-inclusive child bounds,
so a physical-side topology needs no new schema field. The public split-run
case now writes a four-level x-upper tree, stops after the selected-root
checkpoint commit, and reconstructs the same topology in a separate process.

Sparse restart recomputes owners for the new communicator and transfers each
boundary-touching child directly from the I/O root to its selected owner. The
stored geometry and bounds remain authoritative while rank count and ownership
weight stay continuation controls. Every restarted level must still reach the
exact physical side before field parity is accepted.

## Boundary-touching recursive transport (`0.193.0`)

The public fresh and split-run boundary trees now enable Fourier conduction.
Each recursive transport node uses its current fine boundary state on the
x-upper physical side, coarse-time transport context on the remaining sides,
`r^2` child subcycling, time-integrated diffusive registers, physical-side
register omission, reflux, and deepest-to-root average-down.

The sparse path executes the same transport stages on selected owners. The
boundary-touching checkpoint fingerprint records the active transport switch,
thermal conduction, and transport CFL, so continuation cannot silently change
that physics while rank count and ownership weight remain free to change.

## Boundary-touching mixture transport (`0.194.0`)

The same public fresh and restarted trees now activate viscosity,
mixture-averaged species diffusion, and barodiffusion in addition to Fourier
conduction. The recursive transaction therefore carries momentum diffusion,
zero-net-mass correction velocity, species enthalpy flux, and their
time-integrated diffusive registers through every populated level.

The x-upper physical side remains an outflow face: its exterior transport
state comes from the current fine boundary cell and it contributes no
coarse/fine register flux. Every other child side retains coarse-time context,
reflux, and deepest-to-root average-down. Schema 4 already fingerprints all
four transport controls, so this qualification does not require a checkpoint
format change.

## Reacting boundary-touching full physics (`0.195.0`)

The public x-upper trees now enable elementary chemistry together with every
qualified transport term. Each accepted root interval owns one private
candidate across reaction half-step, transport half-step, recursively
subcycled hydro, the second transport half-step, and the second reaction
half-step. No prefix publishes if a later node or stage rejects.

Chemistry is cell-local on the owner of each sparse node, so the physical-side
topology adds no chemistry communication. It does add chemistry state changes
to the coarse/fine transport and hydro contexts and to every average-down. The
schema-4 fingerprint already records chemistry activation, model, and solver
tolerances beside transport controls; rank count and ownership weight remain
continuation-only controls.

## Restart-persistent transport limiter history (`0.196.0`)

Patch-tree checkpoint metadata now stores the minimum accepted transport
limiter theta beside time and minimum root timestep. Serial restart reads it
directly into the accumulated diagnostic before advancing another interval.
Sparse restart reads it only on the selected I/O root, broadcasts it with the
other real metadata, and resumes the communicator-wide minimum from that
value.

The formatted envelope advances from schema 1 to 2 without a fingerprint and
the public fingerprinted format advances from schema 4 to 5. Writers validate
that theta is finite and lies in `[0,1]`; a zero-step checkpoint must retain
the neutral value `1`. A failed or incompatible read publishes the same
neutral value with the existing empty-tree rollback state.

## Restart-persistent conservation baseline (`0.197.0`)

Every patch-tree checkpoint now stores one finite conserved-component vector
representing the composite integral at the beginning of the logical run. The
serial driver restores this vector instead of recomputing it from the restart
state. The sparse driver reads it only on the selected I/O root, broadcasts it
before ownership reconstruction, and uses it for the final run-wide
conservation diagnostic after a changed-rank continuation.

The base envelope advances from schema 2 to 3 and the fingerprinted public
format from schema 5 to 6. A caller that does not supply an explicit baseline
receives a checkpoint-local composite baseline for API compatibility. Public
drivers always pass the original run baseline. An invalid size, nonfinite
component, or MPI rank disagreement rejects the write before file replacement;
a failed read leaves the optional baseline unallocated with the existing
empty-tree rollback state.

## Restart-persistent AMR operator counters (`0.198.0`)

Public patch-tree drivers now accumulate chemistry, transport, and hydro
patch advances in three vectors sized to `patch_tree_maximum_levels`. A
regrid updates only the populated prefix, so removing a deep level does not
erase its earlier work and recreating that level resumes the same cumulative
slot.

The checkpoint envelope advances to base schema 4 and fingerprinted schema 7.
Each file records the common vector capacity followed by the three
nonnegative vectors. The capacity must cover every stored topology level and
fit within the restart configuration. Sparse physics produces owner-local
deltas; the application sums them across the communicator before accumulation
and requires every rank to present identical global vectors to the selected
I/O root. Restart broadcasts all counters before repartitioning numerical
fields. Failed reads publish no optional counter arrays.

## Restart-persistent AMR regrid history (`0.199.0`)

The public patch-tree applications now retain two cumulative adaptation
diagnostics: the number of successful scheduled tag/regrid evaluations and the
sum of cells tagged by those evaluations. An evaluation is counted after the
complete regrid transaction succeeds, whether or not its candidate topology
differs from the current tree. Failed planning, geometry construction, or
migration contributes nothing.

Base checkpoint schema 5 and fingerprinted schema 8 store both nonnegative
integers after the operator-counter vectors. The evaluation count must cover
the committed regrid count and cannot exceed one initialization evaluation
plus one evaluation per committed root step. Sparse writers require exact
communicator agreement before gathering fields; restart broadcasts both values
with the clock metadata before redistributing owners.

## Public branching patch-tree lifecycle (`0.200.0`)

The public boundary and restart inputs now initialize two separated reactive
temperature features. One feature is centered on the x-upper physical side;
the other remains in the interior. Zero-gap clustering therefore produces at
least two ordered patches on a populated level while the boundary branch
continues recursively through all four configured levels.

No numerical representation or checkpoint field changes. Existing ordered
parent/child relations already encode branching, and schema 8 serializes the
complete topology. The stronger fresh and restart gates inspect composite
level/patch identities so a one-patch chain cannot satisfy the public case,
then retain exact serial/sparse and changed-rank field comparisons across both
branches.

## Reproducible completion boundary (`0.201.0`)

The repository freezes its default PeleC/PelePhysics/AMReX/SUNDIALS comparison
snapshot in `references/pelec_baseline.json`. CTest checks that project,
runtime, README, validation, preset, and mapping records remain synchronized.

Serial and sparse patch-tree geometry callbacks now receive an explicit
unlimited-polymorphic context. Production configuration is selected from that
context instead of host association. Callback implementations live at module
scope, so GNU Fortran does not emit stack trampolines. A final-binary gate
parses `PT_GNU_STACK` and rejects executable stack permission.

## Three-dimensional directional core (`0.202.0`)

`mesh_3d_mod` builds independent uniform x, y, and z cell-center arrays using
the established one-dimensional mesh contract. The conserved and primitive
state layouts already carry all three momentum or velocity components.

`directional_flux_mod` now maps z-normal states into the qualified x-normal
Riemann solvers and maps the resulting flux back without changing density,
energy, or transverse momentum. `multispecies_flux_mod` applies the same
rotation while leaving derived thermodynamic and species components in place.
`reactive_directional_flux_3d_mod` performs the corresponding rotation around
the general-EOS reacting Riemann solver, preserving every species flux.

These modules establish coordinate consistency only. There is no 3D field
advance, CFL reduction, boundary fill, transport divergence, AMR hierarchy,
restart format, MPI decomposition, or application at this milestone.

## Conservative uniform-grid 3D Euler path (`0.203.0`)

`finite_volume_3d_mod` stores one periodic cell-centered state and one upper
face flux per cell in each coordinate direction. Each lower face is the
wrapped upper face of the preceding cell, so the x/y/z divergence telescopes
to roundoff without ghost-cell ambiguity. The existing x solver and the
qualified y/z rotations remain the only Riemann implementations.

The multidimensional CFL rate is the per-cell sum of
`(|u|+c)/dx + (|v|+c)/dy + (|w|+c)/dz`. Two PCM spatial evaluations are
composed by SSPRK2. Both intermediate stages live in private candidate arrays;
an invalid face solver or nonphysical stage leaves the caller state unchanged.

`pelef3d` reads a strict periodic namelist, initializes a constant-pressure
entropy wave, clips the last step to the requested time, reports every Euler
integral, and writes deterministic x-fastest CSV output. This architecture is
single-level, serial, constant-`gamma`, and inviscid. It does not yet include
3D reacting state evolution, transport, high-order reconstruction, AMR,
restart, MPI, EB, or production output.

## General-EOS multispecies 3D Euler path (`0.204.0`)

`reactive_3d_mod` extends the direct periodic x/y/z divergence to the complete
runtime reactive state. It calls the established general-EOS x flux and the
qualified y/z momentum rotations, so species components are never reordered.
The multidimensional CFL rate uses the recovered frozen-mixture sound speed.

Each SSPRK2 stage is formed in private state and temperature arrays. Every
cell must pass NASA7 internal-energy inversion, positive pressure and sound
speed, and species-density closure before the next flux evaluation or final
publication. A rejected face solver or candidate leaves both caller arrays
unchanged. `reactive_integrals_3d` integrates every Euler and species component
rather than truncating diagnostics to the first five entries.

`pelef_reactive_3d` selects the elementary seven-species or full ten-species
H2/O2 thermodynamic table, initializes a constant-composition,
constant-pressure diagonal entropy wave, clips the final step, and writes
deterministic x-fastest `Y_k` and `rhoY_k` fields. Chemistry and molecular
transport are deliberately absent from this boundary; their split operators
remain the next regular-grid increments.

## Cell-local chemistry splitting in 3D (`0.205.0`)

`advance_reactive_chemistry_3d` applies the established constant-volume cell
reactor independently to every cell of a private 3D candidate. Each z plane
uses the same transactional 2D chemistry kernel, including the elementary
seven-species explicit path and the full ten-species implicit path. No caller
state or temperature is published unless every plane succeeds.

`advance_reactive_strang_3d` owns a second, outer transaction and composes
`R(dt/2)-H(dt)-R(dt/2)`. Thus a hydro rejection after the first chemistry
half-step or a failure in the final chemistry half-step restores the complete
pre-step state and temperature. The chemistry operator preserves cell density,
three momenta, total energy, and elemental composition while changing species
and temperature through the established reactor contract.

`pelef_reactive_3d` now accepts `entropy_wave`, `uniform_reactor`, and
`reactive_hotspot` initial states plus strict chemistry controls. The public
periodic hotspot initializes a constant-pressure Gaussian temperature field,
advances coupled chemistry and hydrodynamics, and retains deterministic dynamic
species CSV output. Molecular transport, high-order reconstruction, AMR,
restart, MPI decomposition, EB, and external ignition validation remain
outside this boundary.

## Regular-grid molecular transport in 3D (`0.206.0`)

`reactive_transport_3d_mod` evaluates one upper diffusive face per cell in x,
y, and z. It reuses the qualified 2D face coefficient and species-flux
kernels, including mixture-averaged diffusion, barodiffusion, correction
velocity, and species enthalpy. A 3D face kernel reconstructs the complete
velocity-gradient tensor: the normal derivative comes from the two face
cells, while tangential derivatives are centered in their periodic planes.
The Newtonian stress uses Stokes' hypothesis and the energy flux combines
stress work, Fourier heat flux, and species enthalpy transport.

Each cell bounds its total outgoing species mass over all six faces. The
minimum adjacent bound scales each stored face species flux before the
conservative divergence. Two transport Euler candidates form SSPRK2, with EOS
temperature recovery after each candidate and publication only after the
second stage is physical. The parabolic selector uses the largest enabled
momentum, thermal, or species diffusivity and the sum of all three inverse
squared spacings.

`advance_reactive_full_3d` owns a private whole-field candidate and composes
optional chemistry and transport around hydro as
`R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)`. A failure at any stage discards the
candidate and leaves both caller arrays unchanged. `pelef_reactive_3d` loads
the transport table matching the selected elementary or full-H2O2 model,
takes the minimum hydro/transport stable interval, and reports the maximum
diffusivity and cumulative minimum species limiter.

This architecture remains serial, periodic, single-level, and uniform-grid.
High-order hydro, nonperiodic transport boundaries, AMR/reflux, restart, MPI,
and EB are separate increments.

## Static two-level reactive 3D AMR hydro (`0.207.0`)

`amr_hierarchy_3d_mod` defines one rectangular coarse-cell patch and its
integer refinement ratio. PCM prolongation maps each parent value to `r^3`
children, restriction returns their arithmetic volume average, and composite
integrals exclude covered coarse cells before adding fine-cell volume.

The coarse SSPRK2 operator now optionally returns the arithmetic mean of the
two stage face fluxes. A provisional coarse step supplies both endpoint states
for fine ghost construction. Each of the `r` fine SSPRK2 substeps evaluates
its first and second stage boundary flux with linearly interpolated coarse
states at the corresponding substep endpoints. Interior fine faces remain
direct neighboring-cell Riemann solves.

Fine interface fluxes are averaged over `r^2` children per coarse face and
over all `r` substeps. Reflux replaces the provisional coarse interface flux
on the surrounding uncovered cells in all six directions. Covered parents are
then overwritten by the restricted fine state, and every coarse temperature
is recovered from the synchronized conserved field. The outer update owns all
candidate arrays and commits neither level if any flux, EOS recovery, reflux,
or restriction operation fails.

`pelef_amr_reactive_3d` exposes this boundary through paired
`&reactive_3d` and `&amr_reactive_3d` namelists. It is serial, periodic,
static, inviscid, PCM, and single-patch. Chemistry, transport, high-order
reconstruction/prolongation, regridding, MPI ownership, physical coarse
boundaries, and EB remain separate.

## Static 3D AMR hierarchy checkpoint/restart (`0.208.0`)

`amr_reactive_3d_checkpoint_mod` serializes one synchronized two-level
hierarchy in a versioned formatted schema. The header fingerprints every NASA7
species name, molecular weight, temperature range and coefficient, along with
the thermodynamic model, Riemann solver, boundary mode, physics switches,
coarse domain, patch bounds, refinement ratio, physical extents, and CFL. The
payload stores both conserved levels, their recovered temperatures, physical
time, cumulative coarse-step count, initial composite integrals, and maximum
reflux correction.

The reader first validates the complete header and metadata, then loads both
levels into private candidates. It independently recovers temperature,
requires the stored temperatures to agree, verifies coarse/fine average-down
synchronization and the terminal record, and commits caller state and history
only after every gate passes. Missing, truncated, trailing, nonphysical, or
incompatible input therefore leaves all caller arrays and scalar history
unchanged. A restarted driver resumes CFL selection from the stored time and
step and retains the original conservation baseline.

Checkpoint emission occurs only after a complete synchronized AMR root step.
The public ratio-two full-H2O2 case stops at step four and resumes to step
eight; its final coarse and fine CSV files are byte-identical to the
uninterrupted run. The format is currently serial and single-patch. It does
not establish distributed checkpoint layout, dynamic topology, chemistry or
transport source history, EB geometry, or cross-version schema conversion.

## Distributed-slab static 3D AMR hydro (`0.209.0`)

`mpi_amr_reactive_3d_mod` assigns contiguous coarse and fine x-plane slabs to
every rank. A valid distribution gives each rank at least one plane. For an
owned x plane, that rank evaluates the x-, y-, and z-normal faces anchored on
the plane and both SSPRK2 cell candidates. Fixed-order collective sums then
reconstruct the complete face or state candidate on every rank. CFL scans are
owner-local and use a communicator maximum.

The hierarchy is intentionally replicated at this milestone. Once the coarse
and fine flux candidates are identical, every rank applies the already
qualified serial six-face reflux and average-down in the same order. This
preserves byte identity with the serial solver while isolating distributed
arithmetic and collective failure semantics from later sparse storage.

Before a collective whose count depends on a field, all ranks compare fixed-
size patch/layout metadata. Solver controls and the complete NASA7 name and
coefficient table must also agree bitwise. Hierarchy broadcast first proves a
common root and common array extents, receives into private candidates, checks
finite metadata and positive temperatures, and only then publishes the
caller state.

`pelef_mpi_amr_reactive_3d` reads and writes checkpoints and CSV files only on
rank zero. A checkpoint contains no ownership map; after one root reads and
validates it, the transactional broadcast can seed a different valid rank
count. This qualifies deterministic 1/2/4/8-rank arithmetic and changed-rank
restart, but not sparse hierarchy storage, distributed checkpoint files,
scalable I/O, AMR chemistry/transport, dynamic topology, high-order hydro, or
3D EB.

## Characteristic PLM on regular and static-AMR 3D hydro (`0.210.0`)

`reactive_3d_mod` recovers one complete primitive and frozen sound-speed field,
then forms MC- or minmod-limited characteristic slopes independently in x, y,
and z. The y- and z-normal primitive differences are rotated into the same
x-normal characteristic basis used by the qualified one-dimensional kernel.
One common scale protects density, pressure, and every species fraction before
the left and right face primitives are converted back through the NASA7 EOS.
The existing directional Riemann solvers consume those states. Two such
spatial evaluations form the same transactional method-of-lines SSPRK2 update
and publish their arithmetic-mean face flux.

The serial static hierarchy applies that operator to the periodic coarse
level. Each fine substep builds two ghost layers around the patch. A ghost
state first interpolates the coarse conserved endpoints to the fine-stage
time, recovers primitive variables, and then applies component-limited linear
coarse slopes at the ghost cell's fine-center offset. A common positivity
scale and species renormalization precede EOS conversion. The resulting
extended field supplies characteristic PLM states for every fine interior and
coarse/fine interface face. Substep-time and face-area averaging, six-face
reflux, and average-down are unchanged, so the high-order path retains the
existing synchronization and composite-conservation transaction.

The application configuration selects `pcm` or `characteristic_plm` and an
`mc` or `minmod` limiter. Checkpoint schema 2 fingerprints both controls, and
the continuous/restarted PLM public outputs are byte-identical. This milestone
is serial, periodic, inviscid, static, strictly interior, and single-patch for
AMR. The replicated-state MPI AMR kernel remains PCM-only. No CTU transverse
prediction, PPM, physical boundary, AMR chemistry/transport, dynamic topology,
sparse storage, scalable I/O, or 3D EB claim is made.

## Axis-plane embedded boundaries in 3D (`0.211.0`)

`eb_geometry_3d_mod` owns the first 3D cut-cell metric contract. Its analytic
builder places one x-, y-, or z-normal plane in a Cartesian domain and keeps
the coordinate-positive side as fluid. For every cell it stores volume
fraction, normalized fluid centroid, cell classification, the three pairs of
face aperture and tangential centroid data, and the EB area, physical
centroid, solid-to-fluid unit normal, and area-integrated normal. A type-bound
validator checks allocation lower bounds as well as shapes so that the
zero-based face indexing is part of the contract.

For each cell, the validator evaluates

```text
[dy dz (A_x+ - A_x-), dx dz (A_y+ - A_y-), dx dy (A_z+ - A_z-)]
```

and requires it to equal the stored EB normal integral. This discrete
geometric conservation law fixes the wall-normal sign and makes constant
pressure balance independently testable. A plane outside or at a domain bound
collapses to a fully regular or covered geometry. An interior grid-face-
aligned plane is rejected because no positive-volume cut cell would own its
wall metric in this first representation.

`eb_reactive_wall_flux_3d_mod` recovers pressure from the complete NASA7
conserved state. The geometry normal points from solid to fluid, so the wall
flux uses its negative as the fluid-control-volume outward normal. The wall
contributes only pressure momentum flux. Cartesian x/y/z flux differences are
weighted by their open areas and combined with this wall contribution before
division by the fluid cell volume; covered-cell right-hand sides remain zero.

`eb_reactive_hydro_3d_mod` builds PCM Rusanov, HLLC, or PeleC-style fluxes on
open faces. At the six physical domain sides it supplies the adjacent interior
cell on both Riemann sides, giving a zero-gradient boundary. A forward-Euler
candidate updates active cells, recovers temperature cell by cell, leaves
covered cells unchanged, and publishes neither output field if any operation
fails. No redistribution is applied, so the smallest cut-cell fluid volume
still controls stability. General geometry, higher-order EB fluxes, chemistry,
transport, AMR, restart, and MPI are intentionally outside this milestone.

## Conservative flux redistribution for planar 3D EB (`0.212.0`)

`eb_reactive_redistribution_3d_mod` accepts the conservative residual already
divided by each cell's fluid volume. For every cut cell with fraction `kappa`,
it forms a neighborhood from the cut cell and all active face neighbors whose
stored aperture is positive. The neighborhood residual is weighted by fluid
volume fraction. The cut candidate is

```text
kappa R_cut + (1 - kappa) R_neighborhood.
```

The difference from the original volume-weighted cut residual is returned to
the same neighbors. Dividing that excess by their summed volume fraction
makes the global `sum(kappa R)` unchanged component by component, including
overlapping neighborhoods. Regular cells retain their original residual
before receiving excess, and covered cells remain zero.

The module also owns a low-level transactional general-EOS update from a
caller-supplied residual. `eb_reactive_hydro_3d_mod` composes the same operator
after PCM face construction and aperture/wall divergence in
`advance_reactive_eb_redistributed_euler_3d`. The original raw routine is kept
for diagnostics and for a direct instability control in regression tests.

The current neighborhood is local, six-connected, and first order. This is a
qualified FluxRedist boundary, not weighted StateRedist parity. It has no
separate stable-timestep API, multi-step driver, general geometry, high-order
cut-face states, coupled chemistry/transport, AMR, restart, or MPI ownership.

## Zeroth-order weighted StateRedist for planar 3D EB (`0.213.0`)

The second stabilization route operates on the provisional conserved state,
not on its residual. Every small planar cut cell follows the sole nonzero
component of its aperture-difference normal to one full regular receiver. A
target volume fraction of `0.5` is the default. A non-axis normal, missing
receiver, nonregular receiver, or unsupported target is a transactional
failure rather than an implicit geometry generalization.

Each active cell belongs to its own neighborhood and may also receive one or
more small-cell neighborhoods. The receiver count `nrs` partitions shared
receiver volume. The small-cell neighbor weight raises its effective
neighborhood volume to the target; subtracting the same partition from the
receiver's self weight makes the gather and scatter complementary. The
zeroth-order neighborhood mean is then scattered back and divided by each
recipient's `nrs`. Consequently `sum(kappa U)` is invariant component by
component, a uniform active state is unchanged, and standalone covered output
is exactly zero.

`advance_reactive_eb_state_redistributed_3d` forms the provisional state from
a caller residual, redistributes it, and publishes only after NASA7 recovery
succeeds for every active cell. The integrated
`advance_reactive_eb_state_redistributed_euler_3d` shares one private PCM
face-flux and EB-divergence builder with the raw and FluxRedist paths, so all
three compare the same conservative residual.

This is an order-zero, exact-axis-plane contract. It is not the mature 2D
second-order implementation or broad AMReX/PeleC StateRedist parity. General
multi-neighbor geometry, higher-order reconstruction, a public timestep and
application, coupled chemistry/transport, AMR, restart, and MPI remain open.

## Public stabilized planar 3D EB hydro (`0.214.0`)

`reactive_eb_cfl_3d_mod` computes one full Cartesian spectral rate from all
active NASA7 cells and deliberately omits a `kappa` penalty. This timestep is
valid only in composition with a stabilization algorithm, so the public
configuration accepts `flux_redist` or `state_redist` and rejects `raw`.
Covered storage does not participate in EOS recovery or the maximum rate.

`simulation_config_reactive_eb_3d_mod` owns a dedicated single-record input
contract. It requires a strictly interior non-face-aligned x/y/z plane whose
cut sheet has a regular receiver, a valid solver and target, positive density
and pressure, and a finite step/time domain. The geometry builder remains the
authority after input validation.

`reactive_eb_3d_driver_mod` initializes the frozen elementary H2/O2/N2
density sheet, computes fluid-volume integrals and active-cell extrema, and
writes the public EB CSV. Covered cells retain a valid background state for
diagnostics but contribute zero to integrals. The CSV exposes both Cartesian
centers and physical fluid centroids so geometry and state can be checked
without internal module access.

`pelef_reactive_eb_3d` repeatedly selects a CFL step, clips it to final time,
advances one transactional forward-Euler candidate, and publishes only a
successful state/temperature pair. Fluid mass, energy, species, and tangential
momenta are invariant gates. The normal momentum is not gated because the
embedded wall supplies pressure impulse; its change is a reported observable.

The public app remains a serial first-order hydro boundary. It has no
chemistry or transport splitting, characteristic reconstruction, general
geometry, hierarchy synchronization, restart, MPI ownership, or scalable
output.

## Optional elementary chemistry in planar 3D EB (`0.215.0`)

`advance_reactive_chemistry_3d` now accepts an optional three-dimensional
active mask. Each z plane is advanced through the established transactional
2D cell reactor, but no plane is published unless every plane succeeds. The
planar EB driver builds that mask from noncovered cells and verifies covered
state and temperature storage bitwise after every candidate stage.

`advance_reactive_eb_strang_3d` composes the masked elementary reaction
operator with the existing stabilized hydro path as
`R(dt/2)-H(dt)-R(dt/2)`. The complete state and temperature remain private
until both chemistry half steps and the selected FluxRedist or StateRedist
hydro step succeed. Disabling chemistry routes the same candidate directly
through hydro and retains the `0.214.0` StateRedist arithmetic exactly.

The public application exposes a default-off chemistry switch and positive
relative and absolute reactor tolerances. Its reacting case uses only the
frozen seven-species elementary H2/O2/N2 mechanism. H/O/N mole totals are
integrated over physical fluid volume, while mass, total energy, and
tangential momenta remain the hydrodynamic invariant gates. Species densities
are intentionally excluded from the per-component gate when reactions are
enabled and are instead constrained by closure and elemental conservation.

The application still chooses one full-grid CFL step from the state before
the first reaction half step and uses forward-Euler PCM hydro. Therefore this
milestone qualifies the frozen low-feedback public case, not a general
reaction-aware timestep controller or second-order coupled integrator.
Configurable/stiff mechanisms, molecular transport, higher-order EB
reconstruction, general geometry, AMR/reflux/regridding, restart, MPI
ownership, and scalable output remain open.

## Molecular transport in planar 3D EB (`0.216.0`)

`eb_reactive_transport_3d_mod` reuses the regular-grid 3D face constitutive
law but owns the EB transaction. It recovers primitives only for active cells,
evaluates transport on open interior Cartesian faces, multiplies divergence by
stored face apertures, and divides by physical fluid volume. Physical-domain
faces and the adiabatic slip/impermeable embedded wall carry exactly zero
transport flux. The initial transport bound is evaluated on an active-state
sanitized full domain, so arbitrary covered storage cannot change the step.

Before divergence, a physical-inventory limiter computes the outgoing amount
of every species from each active cell. One face factor scales all species
fluxes on that face and adjusts only the species-enthalpy portion of total-
energy flux. Each transport interval uses SSPRK2; both Euler stages pass
through the same order-zero StateRedist and NASA7 recovery path. The final
blend is recovered transactionally and covered state and temperature are
restored bitwise.

`advance_reactive_eb_full_3d` composes one private
`R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)` candidate. Transport-disabled calls
delegate to the `0.215.0` reaction/hydro wrapper, preserving its arithmetic.
Enabled transport requires StateRedist, a complete transport table, at least
one active process, and species diffusion whenever barodiffusion is selected.
The public application chooses the minimum active-state hydro and transport
bounds before the first reaction stage and logs the selected operator sequence,
maximum diffusivity, and minimum species-flux factor.

This remains a serial, single-level, PCM, exact-axis-plane boundary. It does
not include an EB `kappa` penalty beyond StateRedist, nonzero outer transport
fluxes, isothermal/no-slip/catalytic embedded walls, general face-centroid
interpolation, FluxRedist transport, AMR/reflux/regridding, restart, MPI, or a
reaction-aware post-stage timestep controller.

## Transactional checkpoint/restart for planar 3D EB (`0.217.0`)

`reactive_eb_3d_checkpoint_mod` owns a versioned formatted process boundary
for the existing serial planar application. The envelope records the complete
NASA7 species, elementary reaction, and gas-transport tables; the immutable
mesh, exact plane, physics switches, solver, redistribution, tolerances, and
initial-condition controls; the EB metrics; state and recovered temperature;
the initial conserved/L1/H-O-N inventories; and the clock plus cumulative
minimum-step and transport diagnostics. Output and checkpoint paths, cadence,
stop policy, and continuation limits are not model fingerprints.

The reader allocates private state, temperature, inventory, and metadata
candidates. It checks the complete input/model record, target geometry and
shapes, finite and physical active cells, NASA7 temperature consistency, and
the terminal marker before one final publication. No early read mutates the
application state, so malformed or incompatible input has an all-or-nothing
contract.

`pelef_reactive_eb_3d` initializes the deterministic problem and its baseline
inventories, optionally replaces them from one validated checkpoint, and then
continues the same committed-step loop. A scheduled checkpoint is written only
after the entire `R-T-H-T-R` transaction and its state/temperature publication.
Restored cumulative diagnostics make the final run summary comparable to an
uninterrupted process, not merely the final field.

This layer is formatted, serial, single-file validation I/O. It does not own
MPI decomposition, AMR topology, reflux registers, regridding history, general
3D EB geometry, scalable output, or checkpoint schema migration.

## Build-time Cantera mechanism ingestion (`0.218.0`)

`tools/ingest_cantera_mechanism.py` uses a pinned Cantera API to select one
explicit YAML phase and convert supported records into the repository's
normalized mechanism bundle. Ingestion is outside the Fortran runtime. It
rejects non-ideal phases and unsupported thermo, transport, pressure-rate, and
reaction-order families before output publication, and verifies reaction
element balance from the imported species compositions.

`tools/generate_elementary_mechanism.py` consumes the normalized bundle and
emits one module containing reaction construction and rate/Jacobian kernels,
NASA7 and primitive transport loaders, source reaction indices, duplicate
flags, source SHA-256, selected phase, and source/runtime Cantera versions.
`h2o2_full_thermo_mod` remains a compatibility facade, while
`transport_database_mod` converts the generated primitive transport arrays
into its established records. Thus thermo, kinetics, and transport share one
ordered source without introducing a module cycle.

The supported contract is the pinned ten-species, 29-reaction `ohmech` phase:
NASA7 at 101325 Pa, ideal-gas Lennard-Jones transport, elementary and
three-body Arrhenius reactions, and one Troe falloff reaction. The generator
canonicalizes reals with twelve digits after the decimal in exponential form.
Runtime YAML parsing, arbitrary mechanism dispatch, more than 32 species,
CHEMKIN, NASA9, PLOG, Chebyshev, SRI, Lindemann-only falloff, and custom orders
are outside this boundary.

## Build-time selected mechanism bundle (`0.219.0`)

`cmake/MechanismBundle.cmake` turns one normalized JSON bundle into an isolated
build-tree Fortran library. CMake reads the bundle's declared module, loader,
thermo, transport, rate-kernel, Jacobian, and symbol names instead of assuming
the H2/O2 identifiers. A custom command runs the deterministic generator, and
the configured `pelef_mechanism_probe` imports those exact interfaces and
executes them at a temperature inside the common NASA7 interval.

When `PELEF_MECHANISM_SOURCE` is supplied, configuration and every source
regeneration verify its basename and SHA-256 against the normalized bundle.
Both paths are CMake configure dependencies, so a normal rebuild refreshes
bundle-declared interface names and source-hash expectations before generation.
The generated source remains in the build tree and is compiled separately from
`pelef_core`, so the default production targets and committed H2/O2 module are
unchanged. Test builds exercise the same function with a two-species,
one-reaction fixture; a clean tests-disabled build separately qualifies the
pinned ten-species H2/O2 bundle and installs the optional probe.

This boundary proves source provenance, code generation, Fortran module
dependency discovery, loader compatibility, and finite kinetics execution for
one bundle selected at configure time. Generator validation also rejects
nonstandard or colliding Fortran identifiers and emits legal 128-character
reaction equations with bounded-width continuation. This does not route that
bundle into the fixed CFD applications, load mechanisms at runtime, expand the
32-species or supported-reaction limits, or qualify the selected chemistry
physically.

## Build-time selected constant-volume reactor (`0.222.0`)

Selecting a normalized bundle now configures both `pelef_mechanism_probe` and
`pelef0d_selected`. The executable imports the bundle-declared loader, thermo
loader, rate kernel, and species/reaction counts through local aliases. Its
generic runtime dependencies live in `pelef_selected_runtime_core`, whose
Fortran module directory is separate from both `pelef_core` and every
generated mechanism. This prevents the committed full-H2/O2 modules from
colliding with a selected bundle that declares those same module names.

The namelist stores at most 32 `(species name, mole fraction)` entries.
Validation requires finite positive integration controls, finite temperature
and pressure, exact nonduplicate names, a positive composition sum, and blank
and zero unused array tails. Names are mapped exactly into bundle order before
mole-to-mass conversion. Portable lexical input/output aliases are rejected
after reducing `.`, repeated separators, and reducible `..` components.
The initial step must lie within the minimum/maximum bounds and the output
interval cannot be below the minimum. The minimum bounds adaptive reductions;
an output/final scheduler fragment below it is legal only when accepted whole.
Startup then validates each loaded reaction, checks
molecular-mass balance, finds the common NASA7 temperature interval, and
rejects an out-of-range initial state before opening the CSV. The generator
also checks element and molecular-mass balance directly from the normalized
thermo composition records.

The CMake build helper accepts only schema-1 bundles with explicit interface
names, an explicit `chemistry_integrator` value of `explicit` or `implicit`,
and complete NASA7 thermo and primitive transport records aligned one-to-one
with species. The generator emits that policy as a public mechanism constant.
The lower-level generator's legacy kinetics-only mode remains available
independently and is rejected here with a configure-time diagnostic.

Time advancement reuses `advance_constant_volume_implicit_adaptive` with
dynamic NASA7 and elementary-reaction arrays. Density and target specific
internal energy remain fixed. Output columns are generated in mechanism order
for every mass fraction and molar production rate; `wdot_*` uses
`kmol/(m^3 s)`. The generated rate kernel supplies the reported rates, while
the implicit integration remains the shared generic elementary RHS/Jacobian
path.

The native 0D boundary is serial, adiabatic, constant-volume integration for
configure-time selected schema-1 bundles containing 2--32 species and the
already supported elementary, three-body, and Troe falloff forms. It does not
itself provide runtime dispatch or CFD/transport coupling; the separate
regular-1D boundary is described below. Neither boundary establishes physical
validation of arbitrary chemistry.
Symlink, hardlink, mount, case-folding, and
absolute-versus-current-directory path identity remain outside the portable
lexical path check.

### Optional SUNDIALS CVODE backend

`PELEF_ENABLE_SUNDIALS=ON` adds `pelef_selected_cvode_runtime`, separate from
`pelef_selected_runtime_core`. It requires exactly the pinned SUNDIALS 7.2.0
official static Fortran-module libraries for CVODE, serial N_Vector, dense
SUNMatrix, and the dense linear solver. The installed precision header must
declare double precision only; single- and extended-precision ABIs are
rejected before compilation. Shared-only prefixes are rejected so the
installed selected executable does not depend on an unshipped SUNDIALS runtime
library. Fixed applications and the selected mechanism probe do not link this
adapter. The configured selected reactor links it and compiles the CVODE
dispatch only when the option is enabled; otherwise a `cvode` input fails
before output creation.

The CVODE state stores `N-1` independent mass fractions. The largest initial
mass-fraction species is fixed as dependent, so every callback reconstructs
`Y_dep = 1 - sum(Y_independent)`. The RHS uses `reactor_rhs`, and the dense
Jacobian uses `reactor_reduced_jacobian` generalized to the same dependent
index. That Jacobian is energy-constrained and semi-analytic: the kinetics
composition derivative is analytic and `d(rhs)/dT` is centered finite
difference. Temperature remains algebraic and is recovered from the fixed
specific internal energy; it is not a CVODE state component.

The public `cvode_reactor_context` is opaque and stores only a private slot and
generation. A fixed registry owns up to 64 independent sets of Fortran species
and reaction arrays, work arrays, SUNDIALS context, vectors, matrix, linear
solver, CVODE memory, clock, cumulative budget, and failure state. None of
those noninteroperable objects crosses the C callback boundary. Each slot has
a stable `BIND(C)` token containing only its slot and generation. Its address
is installed with `FCVodeSetUserData`; both callbacks validate the token before
resolving the corresponding registry entry.

Live contexts may be advanced in any sequentially interleaved order. The
registry is explicitly not thread safe or reentrant. CVODE trial states that
violate finite composition or closure return a recoverable callback failure.
Caller mass fractions and temperature are committed only after a successful
requested-time return, composition check, energy recovery, and cumulative
internal-step check. A failed context is quarantined until its own
finalization, while other contexts remain usable. Each output call receives
only that context's unused cumulative step budget, so exhaustion is rejected
before an additional internal step.

Finalization first makes the slot unavailable to callbacks, releases CVODE,
the linear solver, matrix, vectors, SUNContext, and Fortran storage in lifetime
order, then clears the public handle. Repeated finalization of that cleared
handle succeeds. An intrinsic handle copy aliases the same live context rather
than copying resources; after the first valid finalization its generation is
stale, so it cannot release a later context that reuses the slot. Exhausting
all 64 slots returns an error without assigning the requested handle. A
nonzero vendor destroy status is terminal: cleanup continues, the handle is
invalidated, and the first failing SUNDIALS destroy routine is reported;
destructor retry is not supported.

`minimum_time_step` normally maps to CVODE's internal minimum. For an output
or final scheduler fragment shorter than that value, the adapter lowers the
minimum only to the requested fragment. `cvode_max_internal_steps` is separate
from the native backend's `maximum_steps`. The public CSV schema and
`wdot_*` units remain identical across backends.

## Build-time selected regular 1D reacting flow (`0.223.0`)

Selecting a normalized schema-1 bundle also configures
`pelef_reactive_1d_selected`. The configured source aliases the bundle's
declared thermo, reaction, transport, and production-rate interfaces and embeds
the normalized bundle SHA-256 for provenance. It loads and validates those
records before opening output, resolves the namelist composition by exact
species name, and supplies the normalized bundle-order mole fractions to the
regular 1D initializer.

The selected executable links `pelef_selected_reactive_1d_runtime` and
`pelef_selected_runtime_core`, never `pelef_core`. The former compiles private
copies of the regular 1D state layout, reconstruction, configuration,
transport, and evolution modules into its own Fortran module directory. The
latter owns generic thermo/kinetics support and the selected generated module.
This separation is required because a selected bundle may deliberately use
the same module and symbol names as the committed full-H2/O2 bundle.

`gas_transport_mod` is the mechanism-independent type and validation layer for
primitive species transport records. `transport_database_mod` reexports that
type while adding the fixed elementary/full-H2/O2 loaders; the selected path
constructs the same type directly from generated names, geometries, and five
transport coefficients. `initialize_reactive_1d` and `simulate_reactive_1d`
accept an optional bundle-order base composition, preserving the fixed-call
source interface when callers are rebuilt and the argument is absent. As with
any changed Fortran module procedure, old objects must not be linked against
the new module ABI. The shared configuration parser accepts
`chemistry_model = "selected"` only when its caller explicitly enables that
mode, so the ordinary fixed application cannot silently reinterpret the new
input.

The resulting path reuses the existing conservative serial regular-grid 1D
hydrodynamics, native chemistry splitting, and mixture molecular transport.
It deliberately rejects AMR and checkpoint/restart settings. CVODE is not
called by this flow solver; selected 2D/3D/AMR/EB/MPI dispatch, runtime
mechanism loading, concurrent use, performance qualification, detailed fuels,
and physical validation remain separate work.

## Build-time selected regular 3D reacting flow (`0.224.0`)

Selecting a normalized schema-1 bundle also configures
`pelef_reactive_3d_selected`. The selected front end enables the otherwise
rejected `thermo_model = "selected"`, loads the bundle-declared thermo,
reaction, and primitive transport arrays, applies the shared structure,
molecular-mass-balance, transport-order, common-temperature-range, and exact
composition-name gates, and then calls `run_reactive_3d_application`.

The fixed `pelef_reactive_3d` front end now performs only fixed model loading
and fixed composition resolution before calling that same application driver.
The driver owns mesh construction, problem initialization, CFL/transport step
selection, transactional `R-T-H-T-R` advancement, invariant checks,
deterministic CSV output, and diagnostics. This makes fixed/selected parity a
comparison of loader/configuration boundaries rather than two copied time
loops.

`pelef_selected_reactive_3d_runtime` links the selected 2D runtime for the
shared thermochemical, regular-face, boundary, and multidimensional core, then
privately compiles the 3D mesh/configuration, directional flux, SSPRK2
chemistry/transport evolution, problem, CSV, and application-driver modules.
Generic transport consumers import only
`gas_transport_mod`; no fixed database loader or committed generated mechanism
enters this link graph.

The selected path therefore supports the existing periodic regular-grid PCM
or characteristic-PLM 3D hydro, native chemistry, and molecular transport for
2--32 selected species. A short nonuniform full-H2/O2
characteristic-PLM/chemistry case keeps that higher-order selected branch live
and byte-identical to the fixed loader. It does not add nonperiodic boundaries,
AMR, EB, MPI, runtime mechanism loading, CVODE inside CFD, thread safety,
performance or detailed-fuel qualification, or external validation.

## Build-time selected regular 2D reacting flow (`0.225.0`)

Selecting a normalized schema-1 bundle now also configures
`pelef_reactive_2d_selected`. The ordinary parser remains fixed-only unless
the selected front end explicitly enables `chemistry_model = "selected"`.
After the shared startup layer validates the generated species, reactions,
transport order, common NASA7 interval, requested initial/hotspot/wall
temperatures, and exact-name input composition, the front end calls
`run_reactive_2d_application` with bundle-order mole fractions.

The fixed 2D executable uses the same application driver after loading its
elementary or full-H2/O2 tables. That driver owns `simulate_reactive_2d`, CSV
publication, extrema, conservation reporting, and all stdout diagnostics.
The initializer, diagonal composition wave, physical-boundary builder, and
simulation accept one optional bundle-order base composition; existing fixed
callers retain their original behavior when it is absent.

`pelef_selected_reactive_2d_runtime` links the selected 1D runtime and
privately compiles `mesh_mod`, the regular 2D configuration and boundary
modules, CTU/physical-boundary evolution, molecular transport, and the shared
application driver. It never links `pelef_core`, either fixed transport
database loader, or either committed generated mechanism. The selected 3D
runtime links this target, so the selected link graph contains only one copy
of the shared regular 2D modules.

The 2D namelist retains bounded 32-species composition and prescribed-flux
storage. The parser validates the complete finite flux record; once the bundle
is loaded, exact resolution additionally requires every slot beyond the actual
selected species count to be zero. This prevents a value intended for an
absent species from being silently ignored at a physical wall.

The supported boundary is serial, single-level regular-grid 2D with the
existing periodic, inflow/outflow, slip/no-slip, adiabatic/isothermal, and
impermeable or prescribed zero-net-mass species-wall controls. Native
chemistry, PCM or characteristic PLM/PPM, optional CTU, and molecular transport
are shared with the fixed executable. AMR, EB, MPI, runtime mechanism loading,
CVODE inside CFD, thread safety, performance, and detailed-fuel validation are
not added by this milestone.

## Build-time selected serial AMR 1D reacting flow (`0.226.0`)

Selecting a normalized schema-1 bundle now also configures
`pelef_amr_reactive_1d_selected`. The front end opts into the otherwise
rejected selected chemistry model, requires AMR, loads and validates the
bundle, resolves exact-name composition, and checks the complete initial
hotspot or constant-pressure entropy-wave temperature range before calling
`run_amr_reactive_1d_application`.

The fixed AMR executable loads its committed model and calls the same driver.
The driver selects the established two-level, arbitrary-depth, or dynamic
multipatch algorithm from configuration, runs it, writes the deterministic
composite CSV, and reports topology and conservation diagnostics. Optional
bundle-order mole fractions reach the root regular-grid initializer through
each AMR entry point. Omitting the argument preserves existing fixed callers.

`pelef_selected_amr_reactive_1d_runtime` links the selected regular-1D runtime
and privately compiles the generic hierarchy, multipatch, regrid, reacting
AMR, and application modules. It never links `pelef_core`, a fixed transport
loader, or a committed generated mechanism. All generic AMR transport
consumers import the mechanism-independent `gas_transport_mod` type.

Selected regular 1D, AMR 1D, 2D, and 3D application drivers pass the generated
bundle integrator constant into their generic chemistry stages. The generic
runtime validates only `explicit` or `implicit`; it does not infer selected
policy from the number or names of species. Fixed front ends omit this optional
argument and retain the pre-existing elementary/full-H2/O2 behavior.

The supported boundary is serial reactive 1D AMR with the existing subcycling,
coarse temporal/spatial boundary fill, reflux, average-down, solution-driven
regridding, native chemistry, and molecular transport. Selected AMR
checkpoint/restart, EB, MPI, runtime loading, CVODE inside CFD, thread safety,
performance, and detailed-fuel validation are not added by this milestone.

## Build-time selected serial reactive EB 2D flow (`0.227.0`)

Selecting a normalized schema-1 bundle now also configures
`pelef_reactive_eb_2d_selected`. The front end opts into the otherwise
rejected selected chemistry model, loads the generated thermo, reaction, and
transport records, resolves exact-name composition, and validates both the
initial-state temperature extrema and any isothermal embedded-wall temperature
before output creation.

Fixed and selected front ends call `run_reactive_eb_2d_application`. This
mechanism-independent driver owns geometry creation, configured outer and
embedded boundary assembly, the single-level time loop, deterministic EB CSV,
volume-weighted invariants, and diagnostics. Optional bundle-order mole
fractions reach both regular-grid initialization and physical-boundary state
construction. The optional generated integrator policy reaches both masked
chemistry half-steps; covered storage remains unchanged and a rejected
chemistry policy leaves the complete caller output unchanged.

`pelef_selected_reactive_eb_2d_runtime` links the selected regular-2D runtime
and privately compiles only the generic 2D EB geometry, reconstruction,
hydro, transport, and application modules. It never links `pelef_core` or a
committed generated mechanism. The public selected target is generated from
the bundle-declared module and loader names and is installed with the other
selected applications.

The supported boundary is serial, single-level reactive 2D EB with the
existing plane/circle geometry, PCM or characteristic PLM, FluxRedist or
zeroth-/second-order StateRedist controls, native chemistry, molecular
transport, and adiabatic/isothermal slip/no-slip embedded walls. Selected EB
outer-domain faces remain outflow-only. Selected EB AMR/3D,
checkpoint/restart, MPI, runtime loading, CVODE inside CFD, general physical
outer boundaries, thread safety, performance, and detailed-fuel validation
are not added by this milestone.

## Build-time selected regular MPI 1D flow (`0.228.0`)

Selecting a normalized schema-1 bundle in an MPI configuration now also
configures `pelef_mpi_reactive_1d_selected`. Fixed and selected front ends load
and validate their own thermo, reaction, and primitive transport records, then
call `run_mpi_reactive_1d_application`. The selected front end additionally
requires the bundle's common NASA7 interval to contain the deterministic
975--1025 K initialization range.

The shared application driver owns the uneven 19-cell periodic decomposition,
bundle-order composition wave, adaptive `R-T-H-T-R` time loop, collective
physicality and conservation checks, ordered gather, and dynamic-species CSV.
The fixed ten-species initialization and column order remain byte-identical to
the pre-extraction implementation. A separate two-species H2/H profile gives
the independent fixture a nontrivial chemistry response.

`pelef_selected_mpi_reactive_1d_runtime` links the selected regular-1D runtime
and `MPI::MPI_Fortran`, privately compiling only the generic MPI domain,
reactive transport, reactive advance, and application-driver modules. It does
not link `pelef_core`, fixed database loaders, or either committed generated
mechanism. The generated explicit/implicit integrator policy reaches both
distributed chemistry half-steps and every adaptive retry. Its integer policy
code must have identical communicator-wide minimum and maximum values.
Invalid or rank-disagreeing policy is reported on every rank without
publishing a candidate state or temperature.

The supported interface remains deliberately narrow: zero or one output-file
argument and one of two exact ordered initialization profiles: H2/H or the
pinned ten-species full-H2/O2 set. It qualifies native chemistry and molecular
transport for those profiles on 1/2/4 ranks. The MPI wrapper and launcher must
resolve to one installation directory, and every fixed/selected coupled MPI
link is audited after build for unresolved or mixed MPI runtime ABIs. This is
not the general serial 1D namelist interface and does not add selected MPI
AMR/EB, physical boundaries, checkpoint/restart, runtime mechanism loading,
CVODE inside CFD, thread safety, scaling, detailed-fuel validation, or external
PeleC field parity.

## Build-time selected serial reactive EB AMR 2D flow (`0.229.0`)

Selecting a complete normalized schema-1 bundle now also configures
`pelef_reactive_eb_amr_2d_selected`. Fixed and selected front ends load and
validate their own models, then call `run_reactive_eb_amr_2d_application`.
Optional bundle-order mole fractions reach configured coarse, fine, and
physical-boundary initialization. The optional generated integrator policy
reaches both chemistry half-steps at both hierarchy levels.

The shared application driver owns the static hierarchy construction,
subcycled coarse/fine schedule, EB reflux, conservative average-down,
deterministic level output, volume-weighted invariants, and diagnostics.
Fixed callers omit both optional arguments and retain the established
numerical path.

`pelef_selected_reactive_eb_amr_2d_runtime` links the selected single-level EB
runtime and privately compiles the generic AMR hierarchy, flux-register,
regrid, transport, driver, and application modules. It does not link
`pelef_core`, fixed database loaders, or either committed generated mechanism.

The supported boundary is serial, static, and exactly two levels. Selected
dynamic regridding, three-level and multipatch modes, checkpoints, and restart
are rejected before output because the current persistence contract does not
record the selected bundle fingerprint and integrator policy. Input, coarse
output, and fine output must be distinct. Initial, hotspot, isothermal-wall,
and reflected-wall extrema must lie inside the bundle's common NASA7 range.
This milestone does not add selected EB 3D, MPI EB AMR, runtime loading, CFD
CVODE, thread safety, performance, detailed-fuel validation, or external
PeleC field parity.

## Build-time selected serial reactive EB 3D flow (`0.230.0`)

Selecting a complete normalized schema-1 bundle now also configures
`pelef_reactive_eb_3d_selected`. Fixed and selected front ends load and
validate their own models, then call
`run_reactive_eb_3d_application`. Optional bundle-order mole fractions reach
regular and cut-cell initialization, and the generated integrator policy
reaches both masked chemistry half-steps in full-step and Strang advancement.
Fixed callers omit both arguments and retain the established seven-species
path.

The shared application driver owns exact axis-plane EB construction, CFL
selection, optional chemistry and molecular transport, hydro advancement,
StateRedist or FluxRedist dispatch, deterministic CSV publication, invariants,
and H/O/N diagnostics. The selected front end additionally checks common
NASA7 temperature bounds and the bounded diagnostic species family before
creating output.

`pelef_selected_reactive_eb_3d_runtime` links the selected regular-3D runtime
and privately compiles the generic EB geometry, CFL, wall flux,
redistribution, hydro, transport, driver, checkpoint-support dependency, and
shared application modules. It does not link `pelef_core`, fixed database
loaders, or either committed generated mechanism.

Selected checkpoint/restart is rejected because the current persistence
schema is fixed-mechanism-specific. The supported boundary is serial,
single-level, static exact non-face-aligned x/y/z planar geometry with PCM,
StateRedist, and the previously qualified transport boundary. This milestone
does not add general embedded geometry or element maps, EB AMR/MPI 3D,
runtime loading, CFD CVODE, thread safety, performance, detailed-fuel
validation, or external PeleC field parity.

## Selected planar EB 3D persistence context (`0.231.0`)

The shared planar EB checkpoint adapter now dispatches two exclusive formatted
schemas without changing the fixed numerical path. A fixed application call
omits every selected-context argument and writes the original schema 1 byte
stream. A selected call must provide all of the configure-time bundle SHA-256,
generated `explicit` or `implicit` integrator policy, and normalized
bundle-order mole fractions; it writes schema 2. Partial context is invalid,
and neither reader accepts the other schema.

Schema 2 writes a `SELECTED_CONTEXT` block immediately after the common header.
The block contains the canonical 64-digit lowercase bundle digest, policy,
species count, and one mole-fraction value per generated species-table entry.
The complete species, reaction, transport, immutable configuration, geometry,
diagnostic, state, temperature, and terminal records that follow remain the
same. The reader compares selected context before allocating state candidates,
then retains the existing rule that no caller-owned geometry, field, clock, or
diagnostic is assigned until every record and end-of-file check succeeds.

`run_reactive_eb_3d_application` chooses the schema from its already validated
fixed/selected call contract. The selected front end permits scheduled writes,
intentional stop, and restart and rejects lexical input/checkpoint or
input/restart aliases. A first-step selected checkpoint continues in a
separate process to the byte-exact uninterrupted final CSV and cumulative
diagnostics.

The bundle digest establishes selected-mechanism provenance only. Schema 2
still relies on record, finite-value, physical-state, temperature, and end
marker validation for payload integrity; it has no cryptographic digest over
geometry, diagnostics, or fields. Finite physically consistent tampering,
crash-atomic file replacement, schema conversion, scalable I/O, general EB
geometry, and AMR/MPI EB 3D restart are not claimed.

## Selected static two-level EB AMR 2D persistence (`0.232.0`)

The shared two-level EB AMR driver now dispatches fixed schema 3 and selected
schema 4 from an all-or-none context contract. Fixed callers omit the bundle
SHA-256, generated chemistry-integrator policy, and bundle-order mole
fractions and preserve the existing byte stream. Selected callers provide all
three values. Partial context and cross-schema reads are rejected.

Schema 4 inserts a `SELECTED_CONTEXT` record directly after the common
four-integer header. The canonical lowercase bundle digest, policy, species
count, and normalized composition precede the established species order,
geometry, embedded-wall identity, numerical controls, patch, clock,
coarse/fine fields, and end marker. The reader validates context and body into
private allocatables and publishes geometry, state, temperature, patch, clock,
steps, and diagnostics only after the terminal record succeeds.

`run_reactive_eb_amr_2d_application` explicitly selects this path for the
selected executable and carries detailed read/write failure context to the
front end. The qualified scope remains serial, static, exactly two levels,
and one fine patch. Dynamic regridding, three-level, multipatch, and MPI EB
AMR persistence do not silently enter schema 4. The front end also rejects
lexical aliases among input, coarse output, fine output, checkpoint, and
restart paths before any read or write.

Schema 4 binds provenance but does not authenticate its body. Crash-atomic
replacement, schema migration, scalable I/O, and cryptographic payload
integrity remain outside the contract.

## Selected dynamic two-level EB AMR 2D persistence (`0.233.0`)

The selected two-level front end now admits the shared serial dynamic-regrid
path when the hierarchy remains exactly two levels with one fine patch. The
checkpoint adapter keeps fixed schema 3 and selected static schema 4
byte-for-byte, then selects exclusive schema 5 for a selected call with
`dynamic_regridding` enabled. A selected static reader expects only schema 4;
a selected dynamic reader expects only schema 5.

Schema 5 reuses the all-or-none selected context record, adds a finite
`nvar`-component `DYNAMIC_BASELINE` record for the initial composite
integrals, and retains the transactional body layout. The body fingerprints
the dynamic policy and stores the actual refined bounds, regrid count, clock,
geometry, diagnostics, and both levels. Therefore restart reconstructs the patch that
exists at the checkpoint rather than the patch originally named in the
namelist. All candidates remain private until the end marker succeeds.

The public full-H2/O2 topology moves from coarse bounds `(2:5,2:5)` to
`(5:12,4:12)` before its first-step checkpoint. Independent continuation
reproduces the uninterrupted three-step coarse and fine CSV bytes and the
uninterrupted selected result equals the fixed dynamic reference. Its maximum
composite conservation error also equals the uninterrupted diagnostic because
restart retains the original baseline. Baseline corruption, composition
mismatch, and static/dynamic cross-schema reads publish no restart target.

This is a serial single-child topology contract. Three-level, multipatch,
dynamic-parent, and MPI EB AMR persistence remain excluded, as do schema
migration, crash-atomic replacement, scalable I/O, and payload authentication.

## Selected static three-level EB AMR 2D persistence (`0.234.0`)

The selected front end now admits the established serial static three-level
path when the hierarchy has one nested patch per level. Fixed static calls
retain the static-three-level magic and schema 3 byte-for-byte. Fixed dynamic-
finest calls retain their separate magic and schema 4 byte-for-byte. A
selected static call uses the static magic with schema 4, which cannot collide
with the fixed dynamic format because the magic strings differ.

The selected schema begins with the all-or-none `SELECTED_CONTEXT` record and
a finite `nvar`-component `COMPOSITE_BASELINE` record. The remaining body
retains the established static topology, geometry, wall and numerical
fingerprints, clock, diagnostics, and all three state/temperature fields. The
reader parses every record into private candidates, validates EOS recovery
and the terminal marker, and only then publishes the root, middle, finest,
patch, clock, diagnostics, and baseline targets.

The selected composition reaches all three initializers and boundaries, and
the generated integrator policy reaches both chemistry half-steps at all
three levels. Restart restores the initial baseline before continuing, so the
cumulative composite-conservation diagnostic is identical to uninterrupted
execution. Front-end lexical checks cover input, all three outputs,
checkpoint, and restart paths.

This is a serial static one-patch-per-level contract. Selected dynamic
three-level, multipatch, dynamic-parent, and MPI EB AMR persistence remain
excluded, as do schema migration, crash-atomic replacement, scalable I/O, and
payload authentication.

## Selected dynamic three-level EB AMR 2D persistence (`0.235.0`)

The selected front end now admits the shared serial dynamic three-level path,
including dynamic-parent regridding, when the hierarchy contains exactly one
middle patch and one finest patch. Fixed calls retain the dynamic-three-level
magic and schema 4 byte-for-byte. Selected calls use schema 5 under that same
magic, so the static and dynamic topology classes remain disjoint.

Immediately after the schema-5 header, the writer emits the all-or-none
`SELECTED_CONTEXT` record and a finite `nvar`-component
`COMPOSITE_BASELINE` record. The remaining bytes follow the established fixed
dynamic body: species order, geometry, wall and numerical policies, dynamic
controls, actual middle and finest patches, clock and regrid diagnostics, all
three state/temperature fields, and the terminal marker. The reader derives
the required schema from the fixed/selected call contract and keeps every
candidate private until context, baseline, topology, fields, recovered
temperatures, terminal marker, and end-of-stream all validate.

The public full-H2/O2 lifecycle moves both nested patches before its first
checkpoint: middle root-cell bounds `(2:11,2:11)` become `(5:10,3:10)`, while
finest middle-cell bounds `(6:9,6:9)` become `(3:10,3:14)`. Restart restores
those committed patches and the original composite baseline. Independent
continuation therefore reproduces the uninterrupted root, middle, and finest
files and its cumulative conservation diagnostic exactly; the uninterrupted
selected result also equals the fixed dynamic reference.

This is a serial exactly-three-level, one-patch-per-level topology contract.
Selected multipatch and MPI EB AMR persistence remain excluded, as do schema
migration, crash-atomic replacement, scalable I/O, and payload authentication.

## Selected dynamic multipatch EB AMR 2D persistence (`0.236.0`)

The selected front end now admits the shared serial dynamic two-level patch-
set path. Fixed calls retain the patch-set magic and schema 3 byte-for-byte.
Selected calls use schema 4 under that same magic. The reader derives the
required schema solely from the fixed/selected call contract, so neither path
can interpret the other format.

Immediately after the schema-4 header, the writer emits the all-or-none
`SELECTED_CONTEXT` record and a finite `nvar`-component
`COMPOSITE_BASELINE` record. The remaining bytes follow the established fixed
patch-set body: species order, geometry, wall and numerical policies, dynamic
controls, complete committed child-patch list, clock and regrid diagnostics,
root state and temperature, every child state and temperature, and the
terminal marker. The reader keeps every candidate private until context,
baseline, topology, fields, recovered temperatures, terminal marker, and
end-of-stream all validate.

The public full-H2/O2 lifecycle advances a 14-by-14 root with two disjoint
10-by-10 child patches. Restart restores both committed patches and the
original composite baseline. Independent continuation therefore reproduces
the uninterrupted root and both child files exactly, and the uninterrupted
selected result also equals the fixed dynamic reference. Bundle-order
composition reaches every initial and boundary state, while the generated
integrator reaches both chemistry half-steps on the root and all children.

The selected front end derives every possible `_patchNNNN` output before
loading the mechanism and rejects lexical aliases with the input, root/base
output, checkpoint, or restart path. This is a serial dynamic exactly-two-
level patch-set contract. Selected MPI EB AMR persistence remains excluded,
as do schema migration, crash-atomic replacement, scalable I/O, and payload
authentication.

## Selected sparse MPI AMR 1D persistence (`0.237.0`)

`mpi_amr_reactive_1d_application_mod` is the mechanism-independent lifecycle
for both fixed and selected sparse MPI AMR 1D front ends. The selected front
end supplies validated generated thermo, kinetics, transport, bundle-order
composition, bundle SHA-256, and integrator policy. Every rank broadcasts and
compares the complete selected context before root initialization, restart, or
output. The sparse advance forwards the selected integrator to both chemistry
half-steps on every owner.

Fixed patch-tree checkpoints retain magic
`PELEF_AMR_PATCH_TREE_REACTIVE_1D_CHECKPOINT` schema 1. Selected calls use
schema 2 and add `SELECTED_CONTEXT` and `COMPOSITE_BASELINE` records around
the established hierarchy body. The reader validates private context,
baseline, hierarchy, finite geometry and clock, nonnegative raw species
densities, EOS-recoverable state, positive finite temperature, terminal
marker, and end-of-stream candidates before publishing. The same physical
state gate precedes selected writes. Failed reads preserve the caller's
solution and baseline; invalid writes leave an existing file unchanged.

The payload contains no rank count or owner map. After read, the application
recomputes distribution for the active communicator, making one-to-two and
one-to-four rank continuation part of the public contract. This is formatted
validation I/O for sparse MPI AMR 1D, not scalable I/O, and selected MPI EB
AMR remains outside the boundary.

## Selected sparse MPI reactive EB patch-tree 2D execution (`0.238.0`)

`mpi_reactive_eb_patch_tree_2d_application_mod` now owns the complete fixed and
selected sparse MPI EB lifecycle: context validation, boundary and root-state
construction, sparse ownership, dynamic regridding, timestep selection,
`R-T-H-T-R` advancement, fixed checkpoint/restart, diagnostics, and composite
CSV output. The two executable front ends only parse configuration, load their
mechanism data, and call this shared driver.

A selected call is all-or-none over bundle SHA-256, bundle-order composition,
and generated chemistry integrator. The driver validates that context and
requires exact communicator-wide agreement before allocating the numerical
root. The single root owner initializes the selected composition, and the
same composition enters boundary construction. The integrator is a trailing
optional argument on the public sparse chemistry, full-physics, and to-time
entry points; full physics forwards it to both chemistry half-steps.

The selected front end and shared driver both reject checkpoint/restart paths,
positive checkpoint cadence, or stop-after-write. This preserves fixed
schema 8 while selected persistence awaits a context-bound format. Fixed
Release checkpoint and stopped-output bytes are frozen against `0.237.0`.
This milestone qualifies selected execution and one-/two-/four-rank parity,
not selected MPI EB persistence, fixed-depth MPI modes, scalable I/O, runtime
mechanism loading, or physical validation.

## Selected sparse MPI reactive EB patch-tree 2D persistence (`0.239.0`)

The shared patch-tree checkpoint adapter now has two exclusive contracts under
the same magic. Fixed calls retain schema 5 or fingerprinted schema 8 exactly.
A call carrying all three selected-context arguments requires fingerprinted
schema 9 and cannot fall back to a fixed reader or writer.

Schema 9 writes species order, `SELECTED_CONTEXT`, bundle SHA-256, generated
integrator, bundle-order composition, the complete fixed numerical
fingerprint, geometry and rank-neutral topology, committed clock,
`COMPOSITE_BASELINE`, original integrals, operator counters, regrid history,
all state/temperature fields, and the terminal marker. Selected read requires
strict end-of-stream after that marker. Context, physical baseline, raw state,
positive temperature, and EOS recovery are validated while the candidate is
private.

The MPI adapter repeats exact context consensus at its public boundary. Root
alone performs formatted I/O; a successful read then broadcasts topology,
derives ownership for the current communicator, and scatters owner-local
fields. Before gather or root read, a collective presence gate requires the
eight selected writer metadata arguments or the selected reader fingerprint
on every rank; rank-local optional-argument returns cannot strand peers in a
collective. Rank count and owner maps never enter the file. Therefore a
one-rank checkpoint can resume on two or four ranks without treating
distribution as physical state.

This remains root-formatted validation I/O. The bundle digest records
mechanism provenance but does not authenticate the payload, and
`status="replace"` is not a crash-atomic commit protocol.

## Rank-local sparse static 3D AMR (`0.240.0`)

`mpi_amr_sparse_reactive_3d_mod` owns the public MPI memory and stepping
boundary for the static periodic two-level 3D hierarchy. The distribution
splits the coarse x extent unevenly and assigns every fine plane to the owner
of its parent coarse plane. Consequently a valid rank always owns coarse
cells but may own a zero-sized fine payload. The hierarchy allocates only the
local conserved and temperature arrays.

The coarse operator propagates two periodic x ghost planes across neighboring
ranks, including through a one-cell local slab, and applies the slab form of
PCM or characteristic PLM. Active fine owners exchange internal x ghosts;
the two exterior sides use coarse start/end halo data interpolated to each
fine SSPRK2 stage and the same limited prolongation as the serial hierarchy.
Fine flux registers are accumulated in deterministic interface order. Reflux,
average-down, synchronized temperature recovery, and the final hierarchy
assignment operate on a private local candidate and publish only after
collective success.

Every public physics entry compares the ownership plan, patch descriptor,
NASA7 species names and coefficients, floating controls, and textual solver/
reconstruction/limiter contract before payload-dependent communication.
Gather additionally compares the full layout before `MPI_Gatherv`. Fine
neighbor ranks are derived from the parent-aligned ownership arrays and are a
distribution invariant.

`pelef_mpi_amr_reactive_3d` keeps global arrays only on the selected root while
initializing, restarting, checkpointing, or writing final CSV. Normal stepping
begins after scatter and deallocation of those arrays. Existing schema 2
contains no communicator size or ownership, so a two-rank checkpoint can
resume at four or eight ranks. This root materialization is compatibility I/O,
not scalable distributed output. The older replicated module remains an
internal arithmetic oracle; dynamic topology, AMR source/diffusion coupling,
physical boundaries, CTU/PPM, and 3D EB AMR are separate responsibilities.

## Fixed chemistry in static 3D AMR (`0.241.0`)

`amr_reactive_3d_mod` now owns a hierarchy source transaction in addition to
hydro. It copies both levels, applies the existing cell-local 3D reactor,
averages fine conserved state over covered coarse cells, recovers the
synchronized coarse temperature, and commits only after every operation is
physical. Its Strang entry composes two such source phases around the existing
hydro transaction. A disabled call dispatches directly to hydro.

`mpi_amr_sparse_reactive_3d_mod` performs the same source operation on local
coarse slabs and any parent-owned fine planes. Empty fine arrays are valid;
their ranks still enter every reaction fingerprint and success collective.
The reaction fingerprint contains equation, type, reversibility, stoichiometry,
all Arrhenius records, third-body efficiencies, and Troe data in source order.
The enclosing split also fingerprints chemistry enablement, tolerances and
integrator policy with the existing solver and reconstruction contract.

Application diagnostics separate the five Euler integrals from reacting
species. H/O/N totals are assembled from composite `rhoY/MW` inventories and
are gated independently. Existing checkpoint schema 2 remains hydro-only;
chemistry-enabled checkpoint or restart input is rejected before stepping
because the schema lacks reaction and source-integrator context.

This boundary is fixed elementary/full-H2/O2, static, strictly interior,
periodic, inviscid, and root-formatted for I/O. Selected mechanisms, transport,
dynamic topology, physical boundaries, 3D EB AMR, scalable I/O, and external
field parity remain separate workstreams.

## Selected mechanisms in static 3D AMR (`0.242.0`)

The serial `amr_reactive_3d_application_mod` and MPI
`mpi_amr_reactive_3d_application_mod` own the numerical lifecycle. Fixed and
generated selected frontends only load their mechanism, resolve the exact
bundle-order composition and output paths, and pass those inputs to the shared
driver. Generated mechanism modules therefore do not enter a fixed runtime
library, and the fixed application does not acquire selected-module linkage.

The sparse MPI driver validates selected provenance before root allocation,
restart, initialization, or file creation. Every rank must agree on optional
context presence, bundle SHA-256, integrator, composition, ordered NASA7
records, ordered reactions including pressure-dependent fields, and ordered
gas-transport records. Real values are compared by representation. The
validator is a public, side-effect-free collective so mismatch tests can
require false on every rank without relying on process abort output.

The lifecycle after validation is unchanged: root initialization, sparse
scatter, rank-local `R-H-R`, gather, invariant checks, and root-formatted CSV.
Ranks without fine planes remain collective participants. H/O/N diagnostics
are enabled only for the qualified H2/O2/N2/AR species-name set; other selected
mechanisms retain physical-state, closure, synchronization, and Euler gates.
Selected transport and persistence are rejected because neither the AMR step
nor schema 2 carries the required selected context.

## Selected static 3D AMR persistence (`0.243.0`)

`amr_reactive_3d_checkpoint_mod` now selects one of two exclusive formats from
the checkpoint call contract. Fixed calls retain schema 2 exactly. Selected
calls require bundle SHA-256, generated integrator, bundle-order composition,
and ordered reactions and write schema 3. Partial selected context returns
before opening the target.

Schema 3 places `SELECTED_CONTEXT` before the established species and state
body. It serializes every reaction equation, type, reversibility flag,
stoichiometric vector, forward/low/high Arrhenius record, Troe record, and
third-body-efficiency vector, followed by chemistry tolerances. The existing
body then binds NASA7, solver/reconstruction/limiter, mesh and patch, CFL,
physics flags, clock, baseline, reflux history, and both levels. Read validates
all records, EOS-recovered temperatures, synchronization, terminal marker, and
strict end-of-stream in private candidates before publishing.

The sparse MPI lifecycle first checks the complete selected context, including
chemistry tolerances, exactly across the communicator. Root performs formatted
read/write, then scatters the restored hierarchy using ownership derived for
the current communicator. No rank count or owner map enters the file, so
changed-rank continuation is part of the contract. This remains root-formatted
validation I/O; transport, payload authentication, crash-atomic replacement,
and scalable I/O are not claimed.

## Static 3D AMR molecular transport (`0.244.0`)

`amr_reactive_transport_3d_mod` composes the established periodic 3D mixture
transport kernel with one static, strictly interior, Cartesian fine patch. A
coarse transport Euler stage advances once while the fine patch advances
`r^2` times. Fine boundary states are linearly interpolated in coarse time and
use limited PLM coarse-to-fine prolongation independently of the selectable
PCM/characteristic-PLM hydro reconstruction. All six fine boundary fluxes are
area averaged and accumulated over the complete parabolic subcycle schedule.

Before reflux replaces a coarse interface flux, a common scalar limiter checks
the uncovered coarse neighbor after removing the old coarse-face contribution.
It bounds species loss with the established 0.9 safety factor, scales the
complete replacement face flux, and applies the equal-and-opposite integrated
correction to the corresponding fine boundary cells. The operation is applied
on x/y/z lower and upper patch faces. Reflux, average-down, and coarse/fine EOS
temperature recovery then complete one private Euler candidate. Two such
candidates form SSPRK2, and the public full step publishes only after the
complete `R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)` transaction succeeds.

`mpi_amr_sparse_reactive_3d_mod` performs the same arithmetic on parent-aligned
rank-local x slabs. Existing coarse and fine halos feed the transport kernel;
time-averaged interface fluxes and per-face limiter values are reduced in a
deterministic layout. Every rank, including ranks with no fine planes,
participates in transport-database, flag, timestep, limiter, reflux, and final
acceptance collectives. Normal stepping never materializes either global
level. Root gather remains confined to compatibility CSV and checkpoint I/O.

The fixed and configure-time selected applications share this lifecycle and
produce exact serial versus one-/two-/four-/eight-rank fields for the qualified
full-H2/O2 case. A transport-disabled call retains the earlier hydro or
`R-H-R` dispatch exactly. Transport-enabled checkpoint and restart remain
rejected because schemas 2 and 3 do not bind the transport database or
parabolic operator history. Dynamic 3D topology, boundary-touching refinement,
physical coarse boundaries, 3D EB AMR, scalable I/O, performance qualification,
and external PeleC field or physical validation remain separate boundaries.

## Static 3D AMR transport checkpoint persistence (`0.245.0`)

The transport-enabled static-AMR checkpoint path uses exclusive schema 4 for
the fixed public `full_h2o2` case with chemistry disabled and schema 5 for the
selected case with `thermo_model='selected'`, chemistry enabled, and transport
enabled. The earlier transport-disabled schema 2/3 contracts remain
byte-compatible and isolated; fixed/selected and transport-enabled/
transport-disabled cross-schema reads are rejected.

The schema-4/5 context includes the complete ordered gas-transport database,
its parameter convention and operator identity, all transport controls, the
`POST_ACCEPTED_COARSE_STEP` phase, and cumulative transport diagnostics. The
checkpoint is written only after the committed
`R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)` hierarchy step. Reads retain all fields,
metadata, recovered temperatures, and diagnostics in private candidates and
publish only after validation, the terminal marker, strict EOF, and close
succeed.

MPI validates selected and transport context on every rank before payload
communication or root I/O. Rank-zero formatted I/O stores no communicator size
or ownership map; successful metadata is broadcast before current-rank
ownership scatter, so exact changed-rank continuation is supported, including
ranks with no fine planes. Root read/write failure is collective. This remains
application-level transactional publication with formatted replacement; it
does not claim crash-atomic replacement, durable commit, payload
authentication, or scalable/distributed I/O. Focused 0.245 evidence and the
known schema-4/5 hashes are recorded in
[`validation/0.245.0.md`](validation/0.245.0.md).
