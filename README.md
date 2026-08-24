# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development`.

## Current capability

The `0.86.0` milestone contains the serial verification suite, seven optional
MPI executables, and runnable serial and sparse-MPI one-dimensional
reactive AMR applications with solution-driven dynamic regridding and
molecular transport. The sparse MPI driver can write an intermediate
patch-tree checkpoint and restart it with a different MPI rank count.

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
Strang sequence applies half reactions only to active cells, performs the EB
hydro transaction, and applies the second half reaction while leaving covered
cells bitwise unchanged. Molecular transport and transverse reconstruction
settings are rejected instead of silently ignored.

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

Unsplit transverse prediction, fourth-order StateRedist slopes,
periodic/ghost-cell neighborhoods, thermal/catalytic wall physics, EB
prolongation and time advancement, dynamic multilevel EB regridding, and MPI
distribution are not yet connected.

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
- direct recursive hydro on sparse AMR payloads with mixed-ratio subcycling,
  replicated flux-register metadata, owner-local reflux/average-down,
  cross-owner PPM face reconciliation, and exact rollback;
- direct recursive molecular transport on sparse AMR payloads with cumulative
  `r²` subcycling, diffusive flux registers, cross-owner shared-face
  reconciliation, and exact rollback;
- a direct sparse `R-T-H-T-R` transaction with owner-only stage execution,
  exact call accounting, serial parity, missing-database rejection, and outer
  rollback after a later-stage failure;
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
