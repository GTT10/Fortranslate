# PeleF architecture

## Executable split

PeleF exposes ten serial verification drivers and six optional MPI drivers
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
  └─ NASA7 general-EOS reactive Euler with PLM/PPM, HLLC, and Strang splitting

pelef_reactive_2d
  └─ NASA7 reactive Euler, physical boundaries, transport, and CTU correction

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
general-EOS Rusanov or HLLC flux
        ↓
conservative finite-volume update

optional transport branch
  ├─ face-centered viscous / conductive / species fluxes
  ├─ explicit SSPRK2 diffusion update
  └─ parabolic dx^2 / diffusivity timestep gate
```

Advective species face fluxes close to the total mass flux. Diffusive species
fluxes use a correction velocity so their sum is zero to roundoff.

## Reactive two-dimensional CTU path

The reactive 2D state is stored as `state(variable,nx,ny)` with a synchronized
`temperature(nx,ny)` field. The normal predictor is selected independently from
the Riemann solver and supports PCM, frozen-composition characteristic PLM, or
time-traced characteristic PPM. For the y direction, momentum and primitive
velocity components are rotated into the x-normal ordering, evaluated by the
same predictor and HLLC/Rusanov kernels, and rotated back.

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
provisional x/y HLLC fluxes
        ↓
conservative transverse half-step correction
  U_face* = U_face - dt/(2 d_t) (F_t,hi - F_t,lo)
        ↓
EOS/positivity bisection on the complete conserved face state
        ↓
final directional HLLC fluxes
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
8. Rusanov remains the robustness baseline. HLLC is the verified
   contact-resolving general-EOS intermediate; it is not labeled as PeleC
   Riemann parity.
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
mirrored velocity/temperature transport gradients, and zero species flux.

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
topology and physical root extent are identical. A deterministic greedy
cell-work schedule then gives each patch exactly one owner; ties go to the
lowest rank. Generic patch fields remain allocated on every rank in this
bridge, but only the owner is authoritative. Collective broadcasts refresh
the replicas, while adjacent sibling faces broadcast the owner's boundary
cells into explicit left/right halo objects for up to four stencil layers.

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
