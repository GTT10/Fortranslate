# Parity and verification strategy

## Verification levels

PeleF uses five gates:

1. unit verification of algebraic kernels;
2. analytical or manufactured-solution verification;
3. parity against a pinned external implementation when available;
4. conservation and element-balance verification;
5. deterministic application-level output checks.

Visual agreement is supplementary and never the sole acceptance criterion.

## Hydro gates

The existing suite retains independent checks for:

- EOS and primitive/conserved conversion;
- Rusanov and qualified PeleC-style Riemann fluxes;
- componentwise and characteristic reconstruction;
- order-2 and order-4 slopes;
- shock flattening;
- smooth entropy-wave convergence;
- Sod exact-solution error;
- Shu-Osher oscillation retention and field signatures;
- planar Sedov-type positivity, symmetry, conservation, and shock location;
- 2D directional flux rotation, CTU transverse corrections, dimensional reduction, and isentropic-vortex convergence.

Higher-order options never replace lower-order baselines in CI.

## Multispecies gates

The passive-species milestone is accepted only when MultiSpecSod reproduces the existing Sod hydrodynamics, each species mass is conserved, `sum_k rho*Y_k` follows `rho`, and 1D/2D smooth species waves converge at approximately second order. A y-uniform 2D multispecies update must reduce to the verified 1D update to roundoff.

## Thermodynamics gates

NASA7 tests pin mass-specific values on both coefficient intervals and verify

```text
cp - cv = R_k
h - u   = R_k T.
```

A fixed O2/N2 mass mixture pins molecular weight, gas constant, `cp`, `cv`, `gamma`, enthalpy, internal energy, pressure, density, and frozen sound speed. Internal energies generated at 300, 1200, and 2500 K must invert back to their source temperatures. Invalid composition and out-of-range energy must fail.

Species molecular weights are aligned with the pinned Cantera 3.2 elemental data used by the runtime parity gate.

## Elementary-kinetics unit gates

The general reaction layer verifies:

- Arrhenius evaluation in SI units;
- valid and invalid stoichiometric records;
- reversible equilibrium constants from NASA7 Gibbs functions;
- forward, reverse, and net progress rates;
- production-rate assembly from net stoichiometry;
- total mass and H/O atom conservation of instantaneous source terms;
- conversion from molar production rates to `dY/dt`.

The generated H2/O2 module must regenerate byte-for-byte from `mechanisms/h2o2_elementary.json`.

## H2/O2 structural reactor gate

The application-level checker reads only the emitted CSV and independently enforces:

- strictly increasing output times and the requested final time;
- constant density;
- finite, non-negative mass fractions;
- mass-fraction closure;
- fixed specific internal energy;
- H and O atom inventories;
- unchanged inert N2;
- instantaneous mass and atom conservation of production rates;
- non-trivial temperature and composition evolution.

## Live Cantera parity

CI installs Cantera 3.2.0 and loads `mechanisms/h2o2_elementary_cantera.yaml`, which describes the same seven species and four reversible reactions as the generated Fortran subset.

Two comparisons are intentionally separated.

### Trajectory comparison

Cantera advances an `IdealGasReactor` from the PeleF initial state using tight solver tolerances. At each PeleF output time, temperature, pressure, and all seven species mass fractions are compared.

### Exact-state production-rate comparison

At each PeleF output row, Cantera is reset to the exact PeleF `(T,rho,Y)` state. Its `net_production_rates` are then compared with the generated Fortran kernel. This prevents integration-history differences from being misidentified as kinetic-rate errors.

Current maximum absolute differences are:

```text
temperature              1.6066e-6 K
pressure                 1.3566e-4 Pa
species mass fraction    1.6967e-11
production rate          3.5527e-12 kmol/(m^3 s)
final temperature        3.6921e-9 K
```

The state relative tolerance is `2e-5`, with an absolute floor of `2e-11`. The rate relative tolerance is `2e-8`, with an absolute floor of `5e-12 kmol/(m^3 s)` to handle nearly cancelled net rates. Thresholds may be changed only with an explained units, data, or numerical-method change.

## Reference-data policy

The repository-wide default is the recursive snapshot in
`references/pelec_baseline.json`, rooted at PeleC commit
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`. A comparison may pin a narrower
source file or dataset, but it must not silently follow a moving branch.

Every external comparison must record:

- upstream repository and commit SHA;
- source mechanism files;
- units and any conversion into SI;
- species ordering and molecular weights;
- initial state and reactor model;
- solver tolerances;
- comparison variables and times.

Pinned numerical signatures may be updated only with an explained method or data change. Conservation limits must not be relaxed merely to accept a regression.

## Reactive one-dimensional gates

The composition-dependent flow path is accepted only when all of the following remain active:

- primitive-to-conserved-to-primitive recovery with nonzero three-component velocity;
- equal-state physical/Rusanov/HLLC/PeleC-style flux identity;
- stationary and moving heterogeneous-composition contact preservation;
- positive general-EOS PeleC-style shock interface state;
- exact equality between summed species flux and total mass flux;
- homogeneous hydro update equal to zero;
- homogeneous Strang-split field equal to independent zero-dimensional cell chemistry;
- global mass, three momenta, and total-energy conservation;
- density, pressure, temperature, and species positivity;
- mass-fraction closure;
- smooth entropy-wave convergence above order 1.75 on both refinement intervals;
- smooth H2/N2 composition-wave convergence above order 1.70;
- smooth H2/N2 composition-wave convergence above order 1.70 with the
  PeleC-style solver selected;
- discontinuous material-contact HLLC error below the Rusanov error;
- nonuniform reactive-hotspot generation of finite pressure and velocity responses.

The hotspot also uses a numerical-reference gate. A 128-cell characteristic-PLM result is restricted onto 32- and 64-cell meshes. At both resolutions, characteristic PLM must have less than 75 percent of the corresponding PCM error, and refinement must reduce the PLM error by at least 30 percent.

This reference is a discretization comparison, not an external physical validation. It establishes that the new high-order path improves on its first-order baseline for the same thermodynamics, reaction model, splitting, and boundary conditions.

## Scope of the evidence

The current Cantera gate establishes parity only for four reversible elementary reactions without third-body or falloff effects. The reactive-flow tests establish numerical coupling and reduction properties for that same subset; they do not establish parity for Cantera's complete `h2o2.yaml`, a stiff mechanism, PelePhysics chemistry integration, full PelePhysics transport, or multidimensional PPM/transport parity. The periodic 2D CTU subset is qualified separately below.

## Reactive PPM gates

The PPM path is compared separately on smooth entropy/composition waves, a moving material contact, and a reacting hotspot. The acceptance criteria require monotonic convergence, conservation, positivity, and sharper material-contact resolution; they do not require PPM to beat characteristic PLM on every smooth problem.

## Reactive characteristic-PPM gates

The time-traced `characteristic_ppm` path remains distinct from the
semidiscrete componentwise `ppm` path.  It is accepted only when:

- five-point reconstruction reproduces linear data and full flattening returns
  the cell center;
- parabolic profile integrals match pinned `u-c`, `u`, and `u+c` reference
  values and reject a characteristic Courant number above one;
- the shock detector returns one in smooth/expanding data and zero for the
  canonical strong compression stencil;
- the contact detector is inactive for pressure jumps and active for the
  canonical density-only contact;
- smooth entropy and composition waves exceed second order on both refinement
  intervals;
- the unsteepened characteristic-PPM material-contact error is lower than the
  componentwise-PPM error;
- bounded contact steepening reduces the characteristic-PPM material-contact
  error by at least 40 percent;
- a periodic pressure-ratio-three shock remains positive and conservative,
  stays inside the initial pressure extrema, and changes when flattening is
  enabled;
- the reacting hotspot remains positive, conservative, and generates finite
  pressure and velocity disturbances.

The characteristic-PPM hotspot is not required to beat characteristic PLM or
componentwise PPM.  Its present value is algorithmic parity of the normal PPM
predictor and sharper contact resolution; full general-EOS characteristic
energy tracing remains outside the claim.


## Reactive two-dimensional CTU gates

The periodic general-EOS CTU path is accepted only while all of the following
remain active:

- equal-state y-direction HLLC flux has the correct normal/tangential momentum
  placement;
- the sum of species fluxes equals the directional mass flux;
- x-normal and y-normal 2D steps agree with the corresponding 1D
  characteristic-PLM and characteristic-PPM updates to roundoff;
- the oblique constant-pressure entropy wave converges under refinement for
  both PLM and characteristic PPM;
- an oblique H2/N2 composition wave converges under characteristic PPM while
  preserving pressure, positivity, and composition closure;
- the transverse correction has a measurable signature and does not materially
  degrade the 32-square result;
- characteristic-PPM shock flattening and bounded contact steepening have
  explicit cell-state signatures;
- an oblique pressure-ratio-three shock remains positive, conservative, and
  free of pressure overshoot while flattening changes the resolved state;
- a periodic 2D material contact is sharper with bounded steepening than
  without it;
- all corrected face states and final cells recover positive density, pressure,
  temperature, and closed composition through the NASA7 EOS;
- global mass, both in-plane momenta, out-of-plane momentum, and total energy
  remain conservative;
- a periodic velocity vortex remains positive and preserves its nontrivial
  velocity/pressure signature;
- both PLM and characteristic-PPM reacting hotspots produce OH and H2O, a
  finite pressure/velocity response, and roundoff-scale composition closure.

The oblique waves are exact advected solutions. The Gaussian velocity vortex is
only a structural regression and is not presented as an analytic steady
solution. The current characteristic-PPM claim covers a PeleC-style normal
parabolic predictor in each coordinate direction followed by the already
qualified conservative full-state CTU half-step correction. It does not claim
complete PeleC multidimensional PPM transverse/corner tracing, physical
boundaries, viscosity, heat conduction, or species diffusion.


## Molecular-transport gates

The one-dimensional transport milestone is accepted only while all of the
following remain active:

- the seven transport records match the thermodynamic species order and retain
  the pinned Cantera Lennard--Jones values;
- pure viscosity, binary-diffusion symmetry, and inverse-pressure scaling are
  fixed unit tests;
- representative mixture viscosity, conductivity, and H2/O2/N2 diffusion
  coefficients are pinned;
- a live Cantera 3.2 probe bounds the known difference between the qualified
  dilute-gas subset and Cantera's mixture-averaged transport model;
- equal left/right states produce zero diffusive flux;
- transverse velocity gradients generate the expected Newtonian shear flux;
- a temperature gradient produces the expected Fourier energy flux;
- diffusive species fluxes sum to zero after correction velocity;
- species enthalpy flux is included in total energy;
- a periodic analytical shear wave converges at second order;
- periodic species and thermal waves smooth while conserving species masses and
  total energy;
- the coupled application remains positive, closed, conservative, and exhibits
  deterministic pressure/velocity signatures.

The Cantera gate is a qualification comparison, not an exact-parity claim. The
allowed envelope reflects the deliberately excluded polar and internal-mode
transport corrections.


## 0.17.0 transport parity gates

The two-dimensional transport path must reduce to the qualified one-dimensional
operator, converge at second order for a periodic shear wave, smooth species and
thermal waves conservatively, and keep the reacting hotspot positive.


## 0.18.0 physical-boundary gates

- the original periodic path retains every 0.17.0 regression;
- slip/no-slip and adiabatic/isothermal ghost states are checked directly;
- solid-wall species flux is exactly zero;
- Couette flow remains linear;
- uniform fixed inflow/outflow remains uniform.

## 0.19.0 chemistry gates

The full mechanism must match its generated source, retain H/O elemental inventories, agree with Cantera in the zero-dimensional reference case, and reduce identically from uniform 1D and 2D CFD states.

## 0.20.0–0.24.0 MPI gates

The distributed one-dimensional milestone is accepted only while all of the
following remain active in both Debug and Release builds:

- uneven 257-cell decomposition and 15-component periodic halo exchange work
  for 1, 2, 4, and 8 ranks;
- global timestep and conservation reductions are communicator-consistent;
- ordered gather reconstructs the same complete field for every rank count;
- passive ten-species Euler transport preserves positivity and species closure;
- eleven independent implicit full-H2/O2 reactors retain energy, closure, and
  nontrivial chemistry response without replicated state;
- distributed general-EOS molecular transport preserves positivity, closure,
  and periodic conservation;
- coupled chemistry, transport, and hydro use a globally synchronized
  accept/reject decision with complete rollback;
- the final coupled field agrees for 1, 2, 4, and 8 ranks within `5e-13`
  relative tolerance;
- the complete serial regression suite still passes in the MPI-enabled build.

## 0.25.0 AMR-foundation gates

The first AMR slice is accepted only while:

- invalid refinement ratios and boundary-touching fine patches are rejected;
- a limited linear profile is prolonged to exact fine-cell averages;
- restricting the prolonged field recovers every covered coarse average;
- average-down changes covered cells and leaves uncovered cells untouched;
- fine substeps sum exactly to the coarse interval;
- coarse and fine interface fluxes accumulate over their respective time steps;
- reflux drives a manufactured composite conservation error to roundoff;
- a used flux register resets so it cannot be applied twice accidentally.

## 0.26.0 AMR-regrid gates

Dynamic regridding is accepted only while:

- a discontinuity tags exactly the adjacent cells of the selected component;
- a flat component produces no tags or fine patch;
- buffering and minimum width produce deterministic strictly interior bounds;
- a physical-boundary tag is rejected instead of silently dropped;
- cells leaving a fine patch are averaged onto the coarse state;
- same-ratio overlap retains the old fine values exactly;
- new fine cells restrict to their source coarse averages;
- patch movement preserves every component's composite integral to roundoff;
- removing and recreating refinement also preserves the composite integral.

## 0.27.0 reactive-AMR gates

The first runnable reactive AMR path is accepted only while:

- the initial reacting hotspot creates a strictly interior fine patch;
- coarse and fine CFL limits produce stable refinement-ratio substeps;
- fine ghost states receive coarse data at each substep time;
- time-integrated coarse and fine interface fluxes are refluxed;
- covered coarse cells equal restricted fine averages after every interval;
- chemistry advances both levels without losing density, momentum, or energy;
- any failed hierarchy operation rolls back both levels transactionally;
- periodic regrid evaluation remains active during the reacting run;
- the final composite state is positive and species mass fractions close;
- composite mass, momentum, and energy remain conserved within `2e-10`;
- CSV cells are coordinate ordered and cover the domain exactly once.

## 0.28.0 reactive-AMR PLM gates

The second-order AMR option is accepted only while:

- limited primitive face states retain positive density and pressure;
- reconstructed species mass fractions are nonnegative and close to one;
- both SSPRK2 stages remain valid on coarse and fine levels;
- their averaged interface flux is the flux accumulated for reflux;
- periodic coarse boundary faces use one identical numerical flux;
- fine substeps receive coarse ghost data at the substep midpoint;
- a smooth moving thermal contact remains dynamically refined;
- PLM has at least 15 percent less composite density error than PCM;
- both PCM and PLM retain composite conservation within `2e-10`.

## 0.29.0 reactive-AMR molecular-transport gates

AMR molecular transport is accepted only while:

- the transport database matches the active thermodynamic species set;
- the coarse timestep includes coarse and `r^2`-scaled fine parabolic limits;
- the fine level completes `r^2` transport substeps per coarse transport
  interval with midpoint-interpolated coarse ghosts;
- coarse/fine gradients use the true adjacent cell-center distance;
- SSPRK2 stage-averaged diffusive fluxes are accumulated for reflux;
- reflux and average-down synchronize every transport half step;
- a conduction-enabled AMR step reduces the temperature span relative to the
  matching inviscid AMR step;
- the dynamic-regrid transport run retains positive thermodynamic states and
  closed, nonnegative species mass fractions;
- every conserved state component and each periodic species mass remain
  conserved within `2e-10`;
- restricted fine cells match covered coarse cells within `5e-13` relative
  tolerance;
- the reacting-hotspot application exercises transport together with
  chemistry, PLM hydro, reflux, and regridding.

## 0.30.0 arbitrary-depth AMR hierarchy gates

The multilevel hierarchy foundation is accepted only while:

- a runtime-sized interface array represents one, two, or arbitrarily many
  levels without a compile-time maximum;
- every adjacent pair has contiguous level numbering, strict nesting, and an
  independently validated refinement ratio;
- cumulative physical bounds, cell counts, and spacings agree across the full
  chain;
- a boundary-touching patch at any depth is rejected;
- recursive prolongation restricts back to every parent cell average;
- cumulative subcycle counts and time steps close exactly to the root interval;
- deepest-to-root average-down synchronizes every covered parent region;
- one independent flux register is owned for each adjacent pair;
- deepest-to-root reflux resets every register and synchronizes each parent;
- a nontrivial four-level, mixed-ratio flux mismatch is corrected to
  roundoff at all three interfaces;
- the arbitrary-depth composite integral is invariant under prolongation,
  restriction, and synchronization.

## 0.31.0 arbitrary-depth reactive AMR gates

The static multilevel reactive engine is accepted only while:

- state and temperature storage is runtime-sized and valid on every level;
- a global root step respects hydro limits from every level after cumulative
  `r` scaling and transport limits after cumulative `r^2` scaling;
- hydro recursively advances each child `r` times with time-interpolated parent
  ghosts and molecular transport advances it `r^2` times;
- every recursive return refluxes the parent and averages the child down before
  the next coarser relation is synchronized;
- chemistry and transport compose as `R-T-H-T-R` across all active levels;
- any failed split operator restores the complete hierarchy transactionally;
- a three-level `24/36/56`-cell hotspot retains the periodic composite integral
  within `3e-10` after PLM hydro and molecular transport;
- every parent covered region equals restriction of its child within `5e-13`;
- all level states retain positive thermodynamics, nonnegative species, and
  mass-fraction closure within `3e-10`;
- a chemistry-enabled three-level split step passes the same synchronization,
  positivity, and closure gates.

## 0.32.0 dynamic multilevel reactive AMR gates

The runnable multilevel regrid path is accepted only while:

- omitting `amr_max_levels` retains the established two-level application;
- a value of three creates levels `0/1/2` from solution-driven tags;
- child-edge tags are excluded so every deeper patch remains strictly nested;
- an unchanged plan retains all existing state without reconstruction;
- a changed plan averages every old level to the root before conservative
  reconstruction and preserves the composite integral within `5e-13`;
- the simulation reevaluates tags at the configured coarse-step interval;
- dynamic three-level periodic hydro conserves mass, momentum, and energy
  within `3e-10`;
- composite CSV rows are spatially ordered, contain all three spacings with
  ratio two, and cover the domain within `3e-13`;
- every output state is finite and retains positive density, pressure, and
  temperature with species closure within `3e-10`.

## 0.33.0 multilevel overlap-transfer gates

Changed multilevel hierarchy transfer is accepted only while:

- old and new physical bounds are intersected independently at each common
  fine level;
- direct transfer occurs only for equal spacing and cell-aligned overlap edges;
- all conserved components and temperature are copied for every aligned cell;
- a spacing mismatch skips direct copy and retains conservative prolongation;
- deepest-to-root average-down follows all copies;
- a forced three-level hierarchy change reports transferred cells;
- the complete deepest overlap is bitwise identical before and after regrid;
- the changed-hierarchy composite integral remains conserved within `5e-13`.

## 0.34.0 multilevel characteristic-PPM gates

The high-order coarse/fine path is accepted only while:

- every level owns four left and right exterior conserved states and
  temperatures;
- physical PPM ghosts preserve periodic wrapping or outflow extrapolation;
- fine PPM ghosts use conservative MC-limited parent subcell values at the
  subcycle midpoint time;
- dynamic PPM patches retain the complete ghost and parent-slope footprint;
- primitive and characteristic PPM advance recursively with SSPRK3;
- flux registers receive the SSPRK3 effective flux used by the state update;
- a three-level characteristic-PPM hotspot conserves the composite state within
  `3e-10`, synchronizes every covered parent region within `5e-13`, retains
  positive thermodynamics and species closure, and completes in Debug/Release;
- the public AMR executable accepts characteristic PPM, produces three ordered
  levels, covers the domain within `3e-13`, and emits only finite positive
  states.

## 0.35.0 multilevel hybrid-WENO5 gates

The WENO slice is accepted only while:

- WENO5-JS and WENO5-Z reproduce constants and linear profiles to roundoff;
- both kernels match fixed nonsymmetric parity points evaluated from PeleC
  `Source/WENO.H` within `3e-14`;
- both schemes advance the fixed three-level reactive hierarchy for two steps;
- composite conservation remains within `3e-10`, covered parent/child states
  synchronize within `5e-13`, temperature remains positive, and species close
  within `3e-10`;
- the public WENO5-Z hotspot case dynamically produces three ordered levels,
  exact composite coverage within `3e-13`, and finite positive states in both
  Debug and Release CI builds.

## 0.36.0 complete hybrid-WENO scheme gates

The remaining upstream schemes are accepted only while:

- WENO7-Z and WENO3-Z reproduce constants and linear profiles to roundoff;
- both kernels match fixed nonsymmetric PeleC `Source/WENO.H` formula points
  within `2e-13` and `3e-14`, respectively;
- schemes 2 and 3 join schemes 0 and 1 in the fixed three-level conservation,
  synchronization, positivity, and species-closure gate;
- the public WENO7-Z and WENO3-Z hotspot cases each dynamically produce three
  ordered levels with exact composite coverage and finite positive states in
  Debug and Release CI builds.

## 0.37.0 physical-boundary AMR gates

Boundary refinement is accepted only while:

- hierarchy, prolongation, restriction, composite integration, and reflux
  accept a patch touching one parent boundary without indexing outside it;
- the physical side receives no reflux and the opposite coarse/fine side
  receives the complete flux-register correction;
- two nested WENO7-Z levels may share the left outflow boundary and advance for
  two steps with composite conservation within `3e-10`;
- every covered parent/child state synchronizes within `5e-13`, thermodynamics
  remain positive, and species close within `3e-10`;
- the public boundary-hotspot case dynamically produces three ordered levels,
  exact coverage within `3e-13`, and finite positive states in Debug/Release.

## 0.38.0 two-level multipatch AMR gates

Multipatch support is accepted only while:

- ordered separated patch sets reject overlap and adjacency before any field
  mutation;
- disconnected tags create multiple plans, while buffered candidates that
  touch or overlap are deterministically coalesced;
- set-wide prolongation, average-down, composite integration, and per-patch
  reflux preserve the parent integral within `5e-12`;
- patch-set movement and repartition preserve the composite integral within
  `5e-12` and retain every aligned same-resolution fine intersection exactly;
- a fixed two-level entropy wave advances two separated WENO7-Z patches for
  two steps with composite conservation within `3e-10`;
- both covered regions synchronize within `5e-13`, thermodynamics remain
  positive, and species close within `3e-10` in Debug and Release CI.

## 0.39.0 multipatch chemistry and transport gates

The reactive multipatch extension is accepted only while:

- requesting transport without a matching database rejects the timestep or
  advance before mutating the solution;
- the root timestep includes every patch hyperbolic limit after `r` scaling
  and every parabolic limit after `r^2` scaling;
- chemistry advances the parent and all patches over equal physical half
  intervals, followed by set-wide average-down;
- transport advances every patch for `r^2` substeps per half interval and
  refluxes one mean-SSPRK2 diffusive flux register per patch;
- a periodic reacting hotspot completes the transactional `R-T-H-T-R` path
  with hydro-quantity conservation within `2e-9`;
- both covered regions synchronize within `8e-13`, thermodynamics remain
  positive, and species remain nonnegative and close within `3e-10`.

## 0.40.0 dynamic two-level multipatch application gates

The tag-driven application integration is accepted only while:

- configuration rejects multipatch mode unless AMR is enabled with exactly
  two levels;
- a deterministic lifecycle starts empty, creates two separated patches,
  moves both patches with nonzero exact overlap transfer, and removes the
  complete set when tags disappear;
- every conserved component retains its composite integral within `8e-12`
  through creation, movement, and removal, with exact regrid counters;
- the public periodic entropy-wave case enables chemistry and molecular
  transport and dynamically creates at least two separated fine segments;
- composite rows are spatially ordered and cover the domain exactly once
  within `3e-13` using level spacings in ratio two;
- density, pressure, and temperature remain positive, while species remain
  nonnegative and close within `3e-10`, in Debug and Release CI.

## 0.41.0 arbitrary-depth multipatch tree gates

The static patch-tree foundation is accepted only while:

- a four-level hierarchy owns `1/2/3/2` patches by level, including a parent
  with no deeper child and two parents that continue refining;
- parent-local children flatten into deterministic level indices and invalid
  parent ownership is rejected;
- physical child bounds remain correct through two ratio-two relations and a
  ratio-three relation;
- recursive conservative prolongation reproduces the root composite integral
  within `5e-12`;
- perturbations on two deepest patches propagate deepest-to-root, synchronize
  every covered parent interval, and retain the composite integral within
  `5e-12` in Debug and Release CI.

## 0.42.0 arbitrary-depth patch-tree synchronization gates

Patch-tree reflux is accepted only while:

- the nested register layout matches every relation, parent, and local child,
  including a parent with an allocated empty register array;
- independent coarse/fine flux mismatches are injected on two deepest patches
  owned by different parents;
- the pre-reflux composite mismatch is nontrivial and deepest-to-root
  synchronization restores the zero reference within `5e-12`;
- every covered interval is synchronized and all registers reset within
  `5e-12` in Debug and Release CI.

## 0.43.0 arbitrary-depth reactive patch-tree hydro gates

Recursive patch-tree hydro is accepted only while:

- a static four-level hierarchy owns `1/2/3/2` reactive patches and contains
  both branching parents and a branch that terminates early;
- one root interval produces exact per-level advance counts of
  `1/4/12/16` under ratio-two subcycling;
- the composite conserved state remains within `3e-10` of its periodic initial
  integral after independent per-child flux-register reflux;
- every parent covered interval matches restriction of its local child within
  `5e-13`;
- every patch retains positive temperature and pressure, nonnegative species,
  and species closure within `3e-10` in Debug and Release CI.

## 0.44.0 reactive patch-tree chemistry gates

Patch-tree chemistry splitting is accepted only while:

- the same four-level branched state advances through a symmetric elementary
  chemistry--recursive-hydro--chemistry interval;
- an otherwise identical hydro-only control proves that chemistry changes at
  least one stored species density by more than `100 epsilon`;
- composite mass, three momentum components, and total energy remain within
  `3e-10` of their pre-step values;
- recursive hydro accounting remains exactly `2/8/24/32` after the second
  accepted root interval;
- every parent-child relation is synchronized within `5e-13`, and all cells
  retain positive temperature and pressure, nonnegative species, and closure
  within `3e-10` in Debug and Release CI.

## 0.45.0 reactive patch-tree transport gates

Patch-tree molecular transport is accepted only while:

- the full reaction--transport--hydro--transport--reaction split advances a
  four-level `1/2/3/2`-patch reactive hotspot;
- two transport half-steps produce exact per-level call counts of
  `2/16/96/256` under ratio-two parabolic subcycling;
- an otherwise identical transport-disabled control differs from the
  transport-enabled state by more than `100 epsilon`;
- composite mass, three momentum components, and total energy remain within
  `2e-9`, and every parent-child relation remains synchronized within `8e-13`;
- temperature and pressure remain positive, species remain nonnegative and
  closed within `3e-10`, and omitting the required transport database rejects
  the step without changing state or counters in Debug and Release CI.

## 0.46.0 runtime patch-tree rebuild gates

Plan-driven patch-tree regridding is accepted only while:

- reapplying an identical four-level plan is a no-op that increments only the
  evaluation counter;
- moving both level-one branches rebuilds all descendants and transfers a
  nonzero number of aligned same-spacing cells by physical coordinate;
- every overlapping deepest state and temperature cell is retained bitwise;
- all composite conserved components remain within `2e-9`, every parent-child
  relation remains synchronized within `8e-13`, and hydro/transport counters
  remain unchanged;
- an invalid parent ownership plan fails without changing state or regrid
  counters in Debug and Release CI.

## 0.47.0 tag-driven patch-tree rebuild gates

Automatic patch-tree planning is accepted only while:

- two separated root features produce two deterministic children owned by the
  root and retain separate branches through three refinement relations;
- per-parent tag clustering reaches the configured four-level limit with
  exact `1/2/2/2` patch counts and deterministic flattened ownership;
- a root-only solution rebuilds transactionally, remains synchronized, and
  preserves every composite conserved component within `3e-10`;
- re-evaluating the unchanged tagged state is a no-op that increments only the
  evaluation counter;
- an invalid tag component fails without changing the solution or counters in
  Debug and Release CI.

## 0.48.0 adjacent patch-tree exchange gates

Independently owned adjacent children are accepted only while:

- patch-tree geometry accepts two touching parent intervals while the older
  separated multipatch API continues to reject them by default;
- each face ghost and all four PPM/WENO ghost layers equal the corresponding
  sibling interior state and temperature exactly before and after advancement;
- a single shared time-integrated hydro flux conserves every composite
  component within `3e-10` under PPM and exact `1/4` level calls;
- the same ownership rule conserves the complete molecular-transport split
  within `2e-9` and produces exact `2/16` transport calls;
- internal fine/fine register sides do not reflux covered parent cells and all
  parent-child relations remain synchronized within `8e-13` in Debug and
  Release CI.

## 0.49.0 MPI AMR patch-distribution gates

The first distributed patch-tree bridge is accepted only while:

- identical valid hierarchy metadata on every rank produces the same owner
  map, and a hierarchy changed on one rank is rejected collectively;
- all nine test patches and all 152 represented cells have exactly one owner;
- communicator sizes 2, 4, and 8 place at least one adjacent sibling face on
  different owners;
- with every non-owner interior poisoned, four halo layers on both sides of
  each adjacent face equal the owner patch's encoded boundary cells exactly;
- broadcasting each owner-authoritative patch reconstructs the complete
  encoded field exactly on every rank;
- the complete gate passes with 1, 2, 4, and 8 ranks in both GNU Fortran Debug
  and Release MPI CI.

## 0.50.0 owner-only MPI AMR chemistry gates

Distributed patch-tree chemistry is accepted only while:

- a four-level `1/2/3/2` tree executes each of its eight patches exactly once
  globally and each rank's call count equals its owner-map patch count;
- the complete distributed reactive state, temperature, and ghost storage
  agree with the serial patch-tree chemistry operator within `5e-13`;
- chemistry changes the reacting state while composite mass, momentum, and
  total energy remain within `3e-10`;
- a deliberately invalid state on a deep patch owner is rejected by every
  rank after earlier owner calls have occurred;
- rejected state, temperature, ghosts, time, and counters return exactly to
  the synchronized pre-call backup;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.51.0 owner-only MPI AMR hydro gates

Distributed patch-tree hydro is accepted only while:

- a four-level `1/2/3/2` PCM tree performs exactly 33 owner patch updates and
  records per-level counts `[1, 4, 12, 16]`;
- each rank's local update count equals the cumulative subcycle weight of its
  owner-map patches, so no patch interval is duplicated or omitted;
- complete distributed state, temperature, ghosts, time, step, and counters
  agree with serial recursive hydro within `5e-13`;
- composite conserved quantities remain within `3e-10` after recursive
  subcycling, shared fluxes, reflux, and average-down;
- a deep invalid owner patch is rejected collectively and the entire solution
  rolls back exactly with zero accepted local calls;
- a two-level PPM tree with six adjacent children crosses at least one MPI
  owner boundary for every multi-rank run, performs exactly 13 updates with
  level counts `[1, 12]`, matches serial within `5e-13`, and conserves within
  `3e-10`;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.52.0 owner-only MPI AMR transport gates

Distributed patch-tree molecular transport is accepted only while:

- a four-level branched tree performs exactly 185 owner patch updates under
  cumulative `r^2` subcycling and records counts `[1, 8, 48, 128]`;
- each rank's count equals the parabolic subcycle weight of its owner-map
  patches, with no duplicated or omitted transport interval;
- complete distributed state, temperature, ghosts, time, step, and counters
  agree with the serial recursive SSPRK2 transport operator within `5e-13`;
- composite mass, momentum, and total energy remain within `2e-9` with
  viscosity, conduction, species diffusion, and barodiffusion enabled;
- a deep invalid owner patch is rejected collectively and the complete tree
  rolls back exactly with zero accepted local calls;
- a two-level tree with six adjacent children crosses MPI owners, performs
  exactly 25 updates with counts `[1, 24]`, matches serial within `5e-13`, and
  conserves within `2e-9` after shared diffusive fluxes and reflux;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.53.0 transactional MPI AMR full-physics gates

The distributed patch-tree `R-T-H-T-R` interval is accepted only while:

- owner synchronization reconstructs root-authoritative time, step, hydro and
  transport counters, and regrid statistics after non-owner poisoning;
- a four-level full-physics step executes exactly 16 chemistry, 33 hydro, and
  370 transport patch calls globally, with each rank matching its owner-map
  cumulative subcycle weights;
- level counters are exactly `[1, 4, 12, 16]` for hydro and
  `[2, 16, 96, 256]` for the two transport half intervals;
- complete distributed fields and bookkeeping match the serial transactional
  patch-tree step within `5e-13` and composite mass, momentum, and total energy
  remain within `2e-9`;
- omitting the required transport database rejects the request before mutation
  and reports zero calls;
- an invalid hydro reconstruction after the first chemistry and transport
  operators rejects the sequence globally, reports zero committed calls, and
  restores the outer backup exactly;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.54.0 sparse MPI AMR storage gates

The rank-local patch container is accepted only while:

- every local sparse patch and cell count equals the distribution metadata;
- global sums contain every patch and cell exactly once;
- the global sum of allocated state, temperature, and ghost values equals one
  replicated tree, independent of MPI rank count;
- gathering into a fully poisoned replica restores fields, ghosts, time,
  steps, level counters, transport counters, and regrid statistics exactly;
- rotating every owner on multi-rank runs reallocates payload only on the new
  owner and preserves the one-copy global value count;
- gathering after migration exactly reconstructs the original replicated
  solution;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.55.0 direct sparse MPI AMR chemistry gates

Sparse patch-tree chemistry is accepted only while:

- every rank advances exactly its owner-map patch count and the global sum is
  exactly one chemistry call per patch;
- deepest-to-root average-down and owner-local temperature recovery reproduce
  the serial four-level branched-tree solution within `5e-13`;
- gathering the sparse result restores state, temperature, normal ghosts,
  PPM-wide ghosts, and bookkeeping without omitted or duplicate payloads;
- a six-adjacent-child PPM tree reproduces serial chemistry within `5e-13`
  while sibling ghost sources cross MPI owners;
- an invalid deepest-level owner state rejects the operation collectively,
  reports zero accepted calls, and restores every local sparse payload exactly;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.56.0 direct sparse MPI AMR hydro gates

Sparse recursive hydro is accepted only while:

- a four-level `1/2/3/2` tree performs exactly 33 owner updates with per-level
  counts `[1, 4, 12, 16]` and each rank matches its cumulative owner weight;
- streamed interval states, effective fluxes, replicated registers, owner-local
  reflux, and average-down reproduce serial hydro within `5e-13`;
- a two-level six-adjacent-child PPM tree performs exactly 13 updates with
  counts `[1, 12]` and matches serial across owner boundaries within `5e-13`;
- gathering after the sparse transaction reconstructs every field, normal and
  wide ghost, time, step, and level counter;
- a negative deepest-level owner state rejects collectively, reports zero
  accepted updates, and restores every sparse payload and counter exactly;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.57.0 direct sparse MPI AMR transport gates

Sparse recursive molecular transport is accepted only while:

- a four-level branched tree performs exactly 185 owner updates under
  cumulative `r²` subcycling with counters `[1, 8, 48, 128]`;
- each rank's call count matches the parabolic weight of its owner patches and
  the gathered state matches serial transport within `5e-13`;
- a two-level six-adjacent-child tree performs exactly 25 updates with counters
  `[1, 24]` and matches serial across shared diffusive owner faces within
  `5e-13`;
- gathering reconstructs state, temperature, ghosts, and transport counters
  after streamed parent data, diffusive registers, reflux, and average-down;
- a negative deepest-level owner state rejects collectively, reports zero
  accepted updates, and restores every sparse payload and counter exactly;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.58.0 sparse MPI AMR full-physics gates

The sparse `R-T-H-T-R` transaction is accepted only while:

- a four-level branched tree performs exactly 16 chemistry, 33 hydro, and 370
  transport owner calls across the communicator;
- per-rank call counts match the owner map and cumulative hyperbolic or
  parabolic subcycle weights;
- counters finish at `[1, 4, 12, 16]` for hydro and `[2, 16, 96, 256]` for
  transport, with exactly one accepted coarse step;
- gathering the sparse result matches the serial full-physics result within
  `5e-13` and preserves composite conserved integrals within `2e-9`;
- a required but missing transport database rejects before mutation and
  reports zero calls;
- a hydro failure after valid chemistry and transport prefixes restores every
  sparse payload, ghost, time, step, and counter exactly and reports zero
  committed calls;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.59.0 topology-changing sparse MPI AMR regrid gates

An explicit-plan sparse topology transition is accepted only while:

- an identical plan leaves fields and ownership unchanged, reports zero
  transferred cells, and increments only the evaluation counter;
- a four-level tree changes patch counts from `[1, 2, 3, 2]` to
  `[1, 3, 3, 2]` and constructs a valid deterministic owner distribution;
- 4/8-rank cases demonstrate that the rebuilt topology changes at least one
  common patch owner, while 1/2-rank cases remain valid without requiring a
  gratuitous reassignment;
- the transfer count and every gathered field, ghost, and counter match the
  serial regrid exactly;
- returned sparse payloads remain globally single-copy and composite
  conserved integrals remain within `2e-9`;
- an invalid parent reference rejects collectively and restores both sparse
  solution and distribution exactly;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.60.0 tag-driven sparse MPI AMR regrid gates

A solution-driven sparse topology transition is accepted only while:

- two separated root momentum features produce deterministic disconnected
  children through four levels with patch counts `[1, 2, 2, 2]`;
- every rank agrees on the positive tagged-cell count, change decision, and
  zero initial fine-overlap count;
- the rebuilt owner payloads, fields, ghosts, counters, and statistics match
  the serial tag-driven regrid exactly;
- the returned hierarchy contains seven globally single-copy sparse patches
  and preserves composite conserved integrals within `2e-9`;
- repeating the same tag decision reports no topology change or transfer and
  increments only evaluation bookkeeping;
- an out-of-range tag component rejects collectively, returns zero counts,
  and restores both sparse solution and owner distribution exactly;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.61.0 point-to-point sparse owner migration gates

Same-hierarchy sparse migration is accepted only while:

- every patch whose owner changes generates exactly one packed send from its
  old owner and one receive on its new owner;
- the 1-rank no-op sends zero messages, while rotated 2/4/8-rank maps report a
  global direct-transfer count equal to the total patch count;
- state, temperature, narrow ghosts, and wide ghosts survive packing and
  unpacking exactly;
- unchanged owners use local assignment without MPI payload traffic;
- the migrated sparse container matches the new owner-map cell and patch
  counts and remains globally single-copy;
- gathering the migrated owners reconstructs the complete original reactive
  tree exactly;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.62.0 point-to-point sparse adjacent halo gates

Sparse adjacent sibling exchange is accepted only while:

- each cross-owner adjacent face generates one bidirectional packed exchange,
  counted once by the left owner;
- the communicator-wide transfer count equals the independently derived
  cross-owner face count, including zero traffic for the 1-rank case;
- PPM sends exactly four state/temperature boundary layers in each direction,
  while narrow reconstruction paths send one layer;
- same-owner siblings copy locally and ranks unrelated to a face allocate no
  send or receive payload;
- the existing six-adjacent-child chemistry, hydro, and molecular-transport
  results retain serial parity within `5e-13` and conservation within `2e-9`;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.63.0 point-to-point sparse child-to-parent gates

Sparse child-interior transfer is accepted only while:

- each child owned by a rank different from its parent generates one direct
  interior-state send to that parent owner during chemistry average-down;
- the communicator-wide transfer count equals the independently derived
  cross-owner child count, including zero traffic for the 1-rank case;
- same-owner pairs copy locally and ranks unrelated to a child/parent pair
  allocate no interior payload;
- hydro and molecular-transport reflux synchronization use the same direct
  transfer helper without changing owner-local subcycle accounting;
- chemistry, hydro, and transport retain serial field parity within `5e-13`,
  composite conservation within `2e-9`, and exact failure rollback;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.64.0 point-to-point sparse parent-state fanout gates

Sparse parent-state fanout is accepted only while:

- each parent owner sends its complete ghost-refresh state once to every
  distinct remote rank that owns at least one child;
- the communicator-wide send count equals an independent owner-map traversal,
  including zero remote recipients for the 1-rank case;
- several children of the same parent and remote owner reuse one received
  parent state;
- same-owner children use the local state and unrelated ranks allocate no
  parent-state payload;
- chemistry, hydro, and transport retain serial field parity within `5e-13`,
  composite conservation within `2e-9`, and exact failure rollback;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.65.0 broadcast-free sparse recursive physics gates

Sparse recursive hydro and molecular transport are accepted only while:

- the sparse physics module contains no `MPI_Bcast` call;
- interval start/end states are packed once per distinct remote child owner
  for every parent invocation;
- each remote child invocation sends one packed left/right boundary-flux pair
  to its parent owner;
- each adjacent shared face sends a correction only for endpoints not owned by
  the parent owner;
- independent hierarchy traversal reproduces all three communicator-wide
  transfer counts under hydro `r` and transport `r^2` subcycle weights;
- level counters equal the serial `[1, 12]` hydro and `[1, 24]` transport
  schedules after one owner-delta reduction per stage;
- adjacent PPM hydro and molecular transport retain serial field parity within
  `5e-13` and composite conservation within their existing tolerances;
- four-level branched parity, exact owner call accounting, and deep failure
  rollback remain unchanged;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.66.0 direct explicit-plan sparse regrid gates

The replica-free explicit-plan topology transition is accepted only while:

- an unchanged plan performs no payload communication and changes only the
  regrid-evaluation counter;
- every new child is conservatively prolongated by its parent owner;
- each cross-owner child receives exactly one packed interior-state payload;
- every retained fine overlap copies state and temperature from its old owner
  directly to its new owner;
- independent old/new hierarchy traversal reproduces the communicator-wide
  prolongation and overlap message counts exactly;
- topology, fields, temperatures, narrow/wide ghosts, counters, and transfer
  accounting match the serial rebuild exactly;
- composite conservation remains within `2e-9`;
- an invalid plan restores the sparse solution and owner map exactly and
  reports zero completed communication;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.67.0 owner-local sparse tag-planning gates

Tag-driven sparse topology construction is accepted only while:

- each candidate parent is tagged exactly once on its owner;
- communicator-wide owner evaluation counts equal the independently derived
  number of evaluated candidate parents;
- only tagged-cell counts and integer child bounds are agreed globally;
- candidate child states remain owner-local and every cross-owner
  prolongation message is counted exactly;
- the final tagged hierarchy matches the serial parent relationships and patch
  extents exactly through four levels;
- the shared direct regrid retains exact state, temperature, narrow/wide ghost,
  counter, and overlap-accounting parity with the serial tagged rebuild;
- installed conserved fields use a deterministic temperature-recovery seed;
- an unchanged tag plan remains a payload-free final commit, while invalid tag
  configuration restores the solution and owner map with zero reported work;
- the sparse MPI module contains neither `MPI_Bcast` nor a materialized-tree
  helper;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.68.0 general-EOS PeleC-style Riemann gates

The reactive PeleC-style acoustic solver is accepted only while:

- equal states reproduce the NASA7 physical flux and expose the same density,
  normal velocity, and pressure at the interface;
- stationary heterogeneous-composition contacts have zero mass and energy
  flux and retain the common pressure;
- moving material contacts reproduce the upwind physical flux, including both
  transverse momenta and every species;
- a nonuniform pressure, velocity, temperature, and composition pair produces
  a finite positive-density/positive-pressure interface state;
- every interface satisfies exact species-flux closure to total mass flux;
- invalid conserved composition is rejected and optional interface outputs are
  reset;
- y-normal equal-state flux agrees with the rotated physical flux;
- the 40/80/160-cell H2/N2 composition wave converges above order 1.70 and the
  160-cell H2 error remains below `4e-6`;
- the complete GNU Fortran Debug and Release suites pass.

## 0.69.0 prescribed species wall-flux gates

The reactive wall transport interface is accepted only while:

- prescribed vectors are finite, match the active mechanism width, occur only
  on slip/no-slip walls, and sum to zero;
- positive input is consistently interpreted from the wall into the gas at
  every lower and upper coordinate face;
- species fluxes close exactly and add their matching NASA7 species-enthalpy
  flux to total energy;
- the positivity limiter scales the entire wall conversion vector and species
  enthalpy contribution together;
- disabled species transport and unbalanced input are rejected;
- a transient lower-wall conversion changes the two selected inventories by
  the integrated boundary flux while preserving total mass and cellwise
  species closure;
- the public 2D application parses and advances the committed namelist case;
- the complete GNU Fortran Debug and Release suites pass.

## 0.70.0 subcycle-weighted MPI AMR ownership gates

MPI AMR work distribution is accepted only while:

- exponent 0 exactly preserves cell-count work, exponent 1 accumulates `r`,
  and exponent 2 accumulates `r^2` across hierarchy depth;
- 64-bit patch work sums exactly to the recorded rank work and global work;
- all ranks reject inconsistent or out-of-range work exponents collectively;
- deterministic least-work assignment does not increase maximum estimated
  parabolic work on the qualification hierarchies;
- the four-level hierarchy strictly reduces maximum estimated work at 2 and 4
  ranks relative to the original cell-only owner map;
- unchanged and changed explicit sparse regrids retain the selected model;
- owner-local tag planning and its intermediate candidate distributions retain
  the same model;
- existing sparse storage, physics, communication-accounting, conservation,
  serial-parity, and rollback gates remain unchanged;
- all gates pass with 1, 2, 4, and 8 ranks in GNU Fortran Debug and Release
  MPI CI.

## 0.71.0 distributed sparse MPI AMR timestep gates

The sparse stability decision is accepted only while:

- every owner evaluates exactly its allocated patch states;
- hyperbolic limits use cumulative `r` scaling and parabolic limits use
  cumulative `r^2` scaling before the global comparison;
- one communicator-wide minimum equals the serial replicated-tree result;
- a root-only hierarchy succeeds with ranks that own no patch payload;
- a nonphysical owner state is rejected by every rank and returns zero `dt`;
- enabled transport without a mechanism-width transport database is rejected
  collectively and returns zero `dt`;
- hydro-only and transport-enabled results pass with 1, 2, 4, and 8 ranks in
  GNU Fortran Debug and Release MPI CI.

## 0.72.0 runnable sparse MPI AMR application gates

The public distributed AMR path is accepted only while:

- a committed namelist reaches a three-level solution-tagged patch tree;
- the selected parabolic work exponent survives initial and periodic regrids;
- stable timesteps are clipped exactly to the configured final time;
- the full `R-T-H-T-R` interval executes on sparse owner payloads;
- final materialization preserves composite mass, momentum, and energy;
- composite CSV cells are strictly ordered and contain each uncovered parent
  or leaf cell exactly once;
- density, pressure, temperature, cell spacing, and species closure are valid;
- CSV headers, row counts, hierarchy levels, coordinates, and all field values
  agree within `5e-13` for 1, 2, 4, and 8 ranks;
- complete GNU Fortran Debug and Release suites pass.

## 0.73.0 rank-independent sparse MPI AMR restart gates

Restart is accepted only while:

- a serial regression round-trips a four-level branching tree, all interior
  conserved states and temperatures, time, steps, regrid accounting, overlap
  accounting, and per-level hydro/transport advances;
- the file identifies its schema, species names, conserved width, base
  geometry, hierarchy depth, refinement plans, and every patch payload;
- a two-rank run writes a three-level intermediate checkpoint and exits cleanly
  before the requested final time;
- four- and eight-rank runs restore that file without using the two-rank owner
  map and finish through the normal sparse physics/regrid loop;
- both restarted composite CSV fields agree with an uninterrupted one-rank run
  within `5e-13` for every header, row, coordinate, and physical value;
- complete GNU Fortran Debug and Release suites pass.

## 0.74.0 embedded-boundary geometry gates

The EB geometry foundation is accepted only while:

- uniformly positive and negative nodal level sets produce exactly regular and
  covered domains with unit and zero cell/face fractions;
- a vertical off-grid plane produces the exact fluid area, one cut-cell column,
  and exact cell and face fractions along that column;
- a diagonal plane produces its analytical fluid area to `3e-13` absolute;
- every geometry validates physical extents, derived spacing, array shapes and
  lower bounds, finite bounded fractions, and classification consistency;
- a circular signed-distance field has nonzero cut cells and its integrated
  area error decreases by more than the required second-order margin from
  `20x20` to `40x40`;
- complete GNU Fortran Debug and Release suites pass.

## 0.75.0 embedded-boundary interface-metric gates

Interface metrics are accepted only while vertical and diagonal planes recover
analytical boundary length, centroid locus, and positive-fluid unit normal to
`3e-13`; a circular interface's perimeter error falls below `0.30` of its
coarse-grid value; every circle normal points inward within `0.02`; and the
refined maximum orientation error is smaller than the coarse error in GNU
Fortran Debug and Release suites.

## 0.76.0 reactive embedded-boundary slip-wall gates

The reactive cut-wall contribution is accepted only while an oblique unit
normal produces `-p*n` momentum flux after general-EOS pressure recovery and
exactly zero mass, z-momentum, energy, and species flux. A 90-degree rotation
must rotate the momentum flux identically, changing velocity at fixed pressure
must not alter it, and vertical and diagonal geometries must integrate to their
analytical pressure forces. Non-unit normals and nonphysical cut-cell states
must be rejected with an all-zero source in GNU Fortran Debug and Release
suites.

## 0.77.0 conservative reactive EB divergence gates

The semidiscrete cut-cell operator is accepted only while it multiplies every
Cartesian face flux by the matching open fraction, adds the integrated
embedded-wall contribution, and divides both by the fluid volume. Uniform
stationary general-EOS pressure must produce zero volume-integrated residual
within `2e-10*p*max(dx,dy)` in regular, vertical-plane, diagonal-plane, and
circular geometries. Covered cells must remain exactly zero, while nonfinite
face fluxes and incorrect face-array extents must fail transactionally in GNU
Fortran Debug and Release suites.

## 0.78.0 conservative EB FluxRedist gates

The first-order small-cell update is accepted only while a cut cell and its
positive-aperture face neighbors form a volume-weighted neighborhood update,
the cut-cell blend removes the inverse-volume stiffness, and the removed
extensive update is distributed with weights whose volume-weighted sum is one.
Every conserved component must retain its domain integral to `4e-13` relative,
a uniform active-cell right-hand side must remain unchanged, covered cells must
remain exactly inert, and redistribution must have one-face support. A
volume-fraction `0.05` reactive cell must remain EOS-valid for a step whose raw
update has negative density; a genuinely nonphysical result must roll back the
entire state and temperature arrays in GNU Fortran Debug and Release suites.

## 0.79.0 weighted StateRedist gates

The zeroth-order StateRedist path is accepted only while its target volume
fraction defaults to `0.5`, merge neighborhoods follow the aperture-normal
ordering, and every recipient is divided by the number of neighborhoods that
contain it. Two diagonal-plane small cells must share one recipient with the
analytical `nrs = 3` weighted result. Both that overlapping case and a full
reactive state must retain every volume-weighted component to `8e-13`
absolute/relative tolerance, and a uniform active state must remain unchanged.
A volume-fraction `0.05` provisional state with negative density must become
EOS-valid, while a stronger nonphysical update, invalid target, or nonfinite
input must fail transactionally in GNU Fortran Debug and Release suites.

## 0.80.0 complete reactive EB hydro gates

The first-order composed update is accepted only if every positive-aperture
face uses the requested reactive Riemann solver, every closed face remains
exactly zero, and nonperiodic domain faces use identical interior states on
both sides. One nonuniform regular-grid face must match the direct HLLC result.
For uniform stationary general-EOS pressure, the complete Riemann-flux,
open-area divergence, integrated-wall-force, weighted-StateRedist, and EOS
pipeline must preserve regular, vertical-plane, diagonal-plane, and circular
geometries within `dt*2e-10*p*max(dx,dy)` per-cell extensive tolerance. Unknown
solvers and nonfinite states must fail with zero fluxes or unchanged state and
temperature in GNU Fortran Debug and Release suites.

## 0.81.0 runnable reactive EB 2D application gates

The public EB path is accepted only while a committed two-namelist input builds
a circular obstacle with regular, cut, and covered cells; the active-cell CFL
matches the analytical uniform-state rate; the clipped time loop reaches the
requested final time; and all conserved components retain their
volume-fraction-weighted integrals. The stationary run must preserve positive
density and pressure, 1000 K temperature, zero velocity, and species closure.
Its CSV must contain 400 ordered grid cells with finite geometry and primitive
fields. Invalid CFL values, geometry without cut cells, and direct requests for
unsupported chemistry must fail without an accepted step in GNU Fortran Debug
and Release suites.

## 0.82.0 active-cell EB chemistry gates

The chemistry-coupled EB path is accepted only while a plane geometry contains
regular, cut, and covered cells; two active-cell reaction half-steps surrounding
the EB hydro update match the regular 2D chemistry application for density,
velocity, pressure, temperature, total energy, and every species field; and the
covered fields remain bitwise identical to a chemistry-disabled EB run. The
reactive result must differ from the inert result, preserve volume-weighted mass
and total energy, and retain all geometry fields. An unknown Riemann solver
after the first reaction half-step must restore the complete input state and
temperature in GNU Fortran Debug and Release suites.

## 0.83.0 active-stencil EB PLM and face-centroid gates

The higher-order EB face path is accepted only while regular and covered faces
have zero centroid offset, a vertical plane gives its analytical partial-face
offset, diagonal-plane partial offsets have the correct sign and magnitude,
and an out-of-range or inconsistent centroid invalidates the geometry. A
synthetic linearly varying face flux must interpolate exactly to each available
partial-face centroid.

On a regular affine moving contact, characteristic PLM must reproduce the
analytical time-centered x- and y-mass fluxes and differ measurably from PCM.
Uniform general-EOS pressure must retain PCM/PLM parity for regular, vertical,
diagonal, and circular geometries through the complete wall-divergence and
StateRedist update. An unknown PLM limiter must return zero face fluxes, and an
invalid hydro selection must retain the caller's state and temperature. The
committed circular-obstacle application must run through the PLM path in GNU
Fortran Debug and Release suites.

## 0.84.0 second-order weighted StateRedist gates

The geometry gate requires zero regular/covered cell-centroid offsets, the
analytical normalized `+0.35` fluid centroid for a vertical-plane cell with
volume fraction `0.30`, bounded finite offsets, and rejection of an invalid
centroid. On overlapping diagonal small-cell neighborhoods, `max_order=2`
must preserve a scalar affine field at every active fluid centroid while the
zeroth-order result remains measurably diffusive. Both smooth and discontinuous
inputs must conserve their volume-weighted scalar exactly; limited recipient
values must remain inside the active input range. Invalid order 1 must return a
zero direct output and must leave a reactive advance state and temperature
unchanged. The committed circular and chemistry applications must run with
`state_redist_max_order=2` in GNU Fortran Debug and Release suites.

## 0.85.0 static two-level EB average-down gates

An aligned fine rectangle must reproduce every covered parent EB volume
fraction by averaging its child fractions. Fine-volume-weighted restriction
must preserve constant states, reproduce an affine scalar evaluated at fluid
centroids on every active parent, leave all outside-patch coarse cells bitwise
unchanged, and make the restricted coarse integral equal the single-count
coarse/fine composite integral to roundoff. A fully covered child block must
exercise the AMReX-compatible first-child fallback.

Uniform multispecies reactive data must restrict without changing conserved
state and must recover temperature on active parents. Covered parents retain
their input state and temperature. A nonphysical active fine block, a nonfinite
input, or inconsistent patch alignment must reject the complete operation
without partial output mutation in GNU Fortran Debug and Release suites.

## 0.86.0 EB flux-register and re-reflux gates

Coarse and fine constant face fluxes representing the same physical interval
must cancel at all four open coarse/fine boundaries when one coarse step is
paired with two fine substeps. A fine-only mismatch crossing a cut interface
must produce a nonzero raw register, scale the cut-cell self correction by its
volume fraction, route a connected recipient below the refined rectangle to
fine children, and preserve the raw fluid-volume-weighted correction in the
post-reflux composite integral.

A successful generic or reactive reflux must reset the register. Nonfinite flux
input must leave accumulated register data unchanged. A uniform reactive state
perturbed proportionally by a small register correction must retain its
temperature after active-cell EOS recovery and preserve covered states exactly.
A correction that makes an active state nonphysical must restore both complete
state arrays, both temperature arrays, and the unconsumed register in GNU
Fortran Debug and Release suites.

## 0.87.0 static two-level reactive EB advance gates

Piecewise-constant prolongation of a uniform multispecies state must inject the
same conserved data into every child, recover active-child temperature, and
return to the original parent state through EB average-down. A proportional
coarse-state change sampled at one-quarter of the coarse interval must produce
the exact interpolated conserved state at every open fine-patch boundary and
must pass EOS recovery; an interpolation time outside `[0,1]` must fail.

A uniform stationary H2/O2/N2 field must survive one characteristic-PLM coarse
step and two fine substeps across a diagonal embedded boundary. Reflux and
average-down must preserve every composite conserved integral, both level
states, active temperatures, and covered data to roundoff. An unknown Riemann
solver must restore both complete state arrays and both temperature arrays in
GNU Fortran Debug and Release suites.

A second hierarchy with a one-percent full-state jump inside the fine patch
must evolve measurably at its coarse/fine boundary while preserving composite
mass, total energy, and every species inventory through fine subcycling and EB
re-reflux. Momentum is not invariant because the embedded wall carries the
pressure reaction force.

## 0.88.0 runnable static reactive EB AMR gates

The public input must construct a 12 by 12 root and an 18 by 18 ratio-two fine
rectangle from the same diagonal level set, with regular, cut, and covered
cells on both levels. A uniform general-EOS mixture must complete the requested
final time in a CFL-valid coarse step, preserve the composite conserved state,
and retain stationary pressure, temperature, velocity, and composition on both
levels.

The installed application must write separate coarse and fine CSV files at the
accepted final time. The checker validates row counts, physical coordinates,
all three EB cell classes, finite positive thermodynamics, stationary fields,
and species closure. Enabling chemistry through the direct static driver must
fail before any accepted step in GNU Fortran Debug and Release suites.

## 0.89.0 solution-driven reactive EB AMR regrid gates

A single-cell temperature hotspot must produce the expected five-point active
tag set while root-boundary and covered cells remain untagged. The buffered
bounding rectangle must remain strictly internal and obey configured minimum
extents. A uniform temperature field must return a valid inactive plan.

Moving a ratio-two patch must first restrict all old fine data, preserve the
complete composite conserved integral, copy every overlapping fine state
exactly, initialize every newly refined active cell by PCM from its synchronized
parent, and leave retired fine regions represented by average-downed root data.
A nonfinite old fine state must reject the transaction without changing root
outputs or publishing a partial new patch.

The public moving-hotspot case must start with a fine patch away from its
temperature maximum, commit at least one regrid, and write a new aligned fine
rectangle that contains the hotspot. Both output levels must retain finite
positive thermodynamics and species closure in GNU Fortran Debug and Release
suites.

## 0.90.0 reactive EB AMR fine-patch lifecycle gates

With initial regrid evaluation disabled, a uniform two-level case must advance
one accepted interval, produce an empty tag plan, conservatively collapse its
fine patch, release all fine arrays and metadata, and finish with a valid
root-only conserved integral. The application must write the synchronized root
CSV and omit the configured fine CSV.

Starting from that root-only state, a new temperature tag must build a valid
strictly internal fine geometry and PCM state without changing the composite
integral. Removing the tag must collapse the re-created patch back to the same
root integral. GNU Fortran Debug and Release suites must also prove that the
inactive driver selects the single-level CFL and hydro path without accessing
unallocated fine storage.

## 0.91.0 reactive EB AMR chemistry gates

An 8 by 8 uniform 1200 K elementary H2/O2 reactor is advanced once through the
regular 2D path and through a ratio-two EB AMR hierarchy whose 10 by 10 fine
patch crosses a plane boundary. Every active coarse and fine EB cell must match
the regular reference density, velocity, pressure, temperature, total energy,
and species mass fractions within `8e-10` scaled tolerance. Both EB levels must
contain regular, cut, and covered cells, the reaction must change temperature,
and active species closure must remain within `5e-13`.

A direct transaction gate lets the first chemistry half-step succeed and then
requests an unknown hydro solver. Failure must leave both full state arrays and
both temperature arrays bitwise equal to their inputs. A second gate advances
chemistry first with an active fine patch and then after conservative patch
collapse through the root-only path; it requires at least two accepted steps
and preserves total mass and total energy within `3e-11` scaled tolerance.
Requesting molecular transport must fail before time, step, or regrid counters
advance in GNU Fortran Debug and Release suites.

## 0.92.0 reactive EB AMR checkpoint/restart gates

Direct round trips must preserve the complete coarse and fine conserved arrays,
actual patch bounds, lifecycle flag, time, accepted-step and regrid counters,
minimum timestep, and base density. Active temperatures recovered from the EOS
must remain within `3e-12` scaled tolerance. A second round trip begins after a
reacting fine patch collapses and must restore a valid root with no allocated
fine state, temperature, geometry, or patch metadata. A file truncated after a
valid magic header must leave all state outputs unallocated and both geometries
invalid.

The public split-run gate starts a uniform 1200 K elementary H2/O2 hierarchy
with one active fine rectangle. After the first reacting interval, the regrid
transaction removes the untagged fine patch, writes a root-only checkpoint, and
stops before final time. Restart must continue from the stored time and counters
through the root-only chemistry/hydro path. Every final CSV value must match an
uninterrupted reference within `3e-10` scaled tolerance, chemistry must progress,
species closure must hold, and no inactive fine CSV may be written.

## 0.93.0 reactive EB AMR multipatch kernel gates

Planner gates require two disconnected tag clusters to produce two ordered
rectangles independent of tag insertion order. Empty tags produce an empty
collection, a configured tag-gap joins nearby components, candidates within
the two-cell EB redistribution safety region coalesce, and any rectangle that
cannot stay strictly internal is rejected.

A ratio-two plane-EB hierarchy uses two separated 10 by 10 fine patches over a
14 by 14 root. Topology gates require set-wide average-down and removal to
preserve the composite integral, movement to retain every aligned old/new fine
cell exactly, unchanged children to remain bitwise identical, newly exposed
cells to equal PCM parent injection, and a NaN candidate to roll back the full
root and patch set.

The hydro gate advances the root once and subcycles both children with distinct
interface fluxes. It requires valid active thermodynamics, exact post-step
average-down synchronization, and composite mass, total-energy, and species
conservation within `5e-11` scaled tolerance. The chemistry gate adds elementary
H2/O2 reaction half-steps on the root and both children, requires a measurable
species change and mass/energy conservation within `5e-10`, and proves final
synchronization. An invalid hydro solver after the first chemistry half-step
must restore every coarse/fine state and temperature bitwise in GNU Fortran
Debug and Release suites.

## 0.94.0 reactive EB AMR multipatch application gates

The public EB AMR driver must convert its configured seed rectangle into a
valid patch set, optionally replace it from disconnected temperature-gradient
tags, and select the minimum root/all-child CFL timestep. Periodic planning
must retain an identical collection, conservatively replace a changed
collection, and remove every child only when the configured empty-tag policy
allows it. Time, accepted-step, minimum-timestep, and regrid accounting may
change only after a complete set-wide physics or topology transaction commits.

An input-driven double-hotspot case over a plane EB must produce two ordered,
separated 10 by 10 ratio-two children over a 14 by 14 root. The public run must
advance chemistry and hydrodynamics, preserve composite mass and total energy,
retain finite positive thermodynamics and species closure, and exercise a
periodic regrid. Output gates require one synchronized root CSV and exactly two
fine CSV files whose deterministic patch suffixes, coordinates, row counts,
temperature features, EB classes, and simulation time match the committed
hierarchy.

Configuration must reject multipatch checkpoint or restart requests before
initialization because formatted checkpoint schema one owns at most one fine
rectangle. The focused public-driver and input/output gates must pass before
the complete GNU Fortran Debug and Release regression suites.

## 0.95.0 reactive EB AMR multipatch checkpoint/restart gates

A direct round trip must preserve the root state, ordered child count, every
child's actual coarse bounds and complete state, time, accepted-step and regrid
counters, minimum timestep, and base density. EOS-recovered active
temperatures must match the written hierarchy within `3e-12` scaled tolerance.
A file truncated after the valid patch-set magic must leave root arrays,
geometry, and the child collection unpublished.

The public split-run gate advances a reacting double-hotspot plane-EB hierarchy
with two separated children. A scheduled checkpoint after the first committed
step and periodic regrid must stop before final time. Restart must reconstruct
both children and continue from the stored counters. Root and both child CSVs
must match an uninterrupted reference in every field within `3e-10` scaled
tolerance, retain finite positive thermodynamics and species closure, preserve
deterministic child order, and reach the requested final time.

The focused direct and public split-run gates must pass before all 175 tests in
GNU Fortran Debug and Release configurations. Existing single-patch checkpoint
round trips and public restart parity remain in the same complete regression.

## 0.96.0 reactive EB AMR physical-boundary patch gates

The public configured-patch case places a ratio-two fine rectangle on the
x-lower outflow boundary of a 12 by 12 plane-EB root. Its fine mesh must begin
exactly at the global lower coordinate, span 12 by 18 cells, and contain
regular, cut, and covered EB classes. The root and fine CSV files must reach
the requested final time with exact mesh-center coordinates, finite positive
thermodynamics, species closure, stationary pressure and temperature, and
negligible velocity drift from a uniform state.

The run exercises the physical-side fine-state exterior closure while retaining
coarse-time interpolation, flux accumulation, reflux, and average-down on the
three coarse/fine sides. The focused nine-test EB AMR application gate must
pass before all 177 tests in GNU Fortran Debug and Release configurations.
Existing planner tests continue to reject physical-boundary tag collections,
making the qualified scope a configured static single patch rather than
dynamic boundary refinement.

## 0.97.0 reactive EB AMR dynamic physical-boundary gates

Unit gates place a temperature jump on a root boundary cell and require a
one-sided tag plus a domain-inclusive single-patch plan. Direct collection
planning must accept a boundary tag. The two-child hydro and Strang gate moves
its first plane-EB child onto the x-lower outflow side, requires composite mass,
energy, and species conservation, exact synchronization, and whole-hierarchy
rollback, then moves that child inward while retaining every aligned overlap
cell exactly and filling only newly exposed cells from the synchronized root.

The public case starts with one hotspot beside x-lower and a second interior
hotspot. Initialization and accepted-step planning must produce two ordered,
separated children, with the first exactly aligned to x-lower. Root and child
outputs must reach the requested time with finite positive thermodynamics,
species closure, ratio-two spacing, aligned inferred bounds, all EB classes,
retained hotspot temperatures, and the established two-cell separation. The
focused eleven-test application gate must pass before all 179 tests in GNU
Fortran Debug and Release configurations.

## 0.98.0 static three-level EB synchronization gates

A root field, a ratio-two middle rectangle, and a ratio-two finest rectangle
nested inside that middle level must produce one composite volume-weighted
integral. After deepest-to-root average-down, integrating the synchronized root
alone must match that composite value while cells outside each child rectangle
remain unchanged. A constant finest field must restrict exactly into its
covered middle region.

The reactive gate applies distinct conservative scalings on all three levels,
requires the same composite conservation for every state component, and
requires finite positive EOS-recovered root and middle temperatures. A
nonfinite finest state or nonpositive finest temperature must reject the whole
operation and return both parent fields unchanged. The focused two-test EB
hierarchy gate must pass before all 180 tests in GNU Fortran Debug and Release
configurations.

## 0.99.0 static three-level reactive EB hydro gates

A ratio-two middle rectangle and a ratio-two finest rectangle produce one root
update, two middle updates, and four finest updates. The finest rectangle is
two middle cells from the middle boundary and its complete interface has unit
open-area fraction. Distinct conservative scalings on the three levels must
evolve through both nested interfaces while composite mass, total energy, and
every species remain conserved. Root and middle must finish synchronized from
the deepest level, and all three temperature fields must remain finite and
positive.

An unknown Riemann solver must leave all three state and temperature fields
unchanged. Mutating one finest boundary face to a fractional EB area must
reject the hierarchy before advancement and retain every input field. The
focused multipatch/multilevel gate must pass before all 180 tests in GNU
Fortran Debug and Release configurations.

## 0.100.0 EB-cut nested-interface conservation gates

The finest rectangle must return to the plane-EB geometry whose boundary has
at least one fractional open-area face. Distinct root, middle, and finest
states must complete the same `1/2/4` recursive hydro schedule while preserving
composite mass, total energy, and every species. All levels must retain finite
positive temperatures and the root and middle must finish deepest-first
synchronized.

The closure must derive its target from the pre-update middle/finest composite
integral and signed middle exterior flux. It may correct only density, total
energy, and species over uncovered active middle cells; EB wall momentum is not
a zero-residual quantity. Solver failure and closure/EOS rejection must retain
all three original state and temperature fields. The focused
multipatch/multilevel gate must pass before all 180 tests in GNU Fortran Debug
and Release configurations.

## 0.101.0 static three-level reactive EB Strang gates

The qualified cut-interface hierarchy receives one chemistry half-step on
each of the root, middle, and finest levels, the established `1/2/4`
recursive hydro schedule, and a second chemistry half-step on every level.
Post-chemistry reactive average-down must make both parent levels identical
to an explicit deepest-first synchronization.

The composite hierarchy must preserve mass and total energy, retain species
sum equal to density, show a nonzero reaction-driven change in at least one
species integral, and keep every temperature finite and positive. Selecting
an unknown Riemann solver after the first chemistry half-step must reject the
whole transaction and return all three state and temperature fields exactly.
The focused multipatch/multilevel gate must pass before all 180 tests in GNU
Fortran Debug and Release configurations.

## 0.102.0 public static three-level EB AMR gates

The public namelist case builds an 8 by 8 root, a ratio-two 12 by 12 middle
rectangle, and a nested ratio-two 16 by 16 finest rectangle whose interface
is crossed by a plane EB. The time loop must choose a positive CFL limit from
all three levels, reach the requested final time, preserve composite mass and
total energy, and emit one CSV per level.

Every CSV must have the expected cell count, finite values, the final time,
regular/cut/covered EB classes, positive thermodynamics, and species closure.
Active root, middle, and finest cells are compared with a uniform regular-grid
reactor reference, which must itself show nonzero chemistry evolution. The
focused public application gate must pass before all 183 tests in GNU Fortran
Debug and Release configurations.

## 0.103.0 static three-level checkpoint/restart gates

An uninterrupted public run and a split run use the same 8 by 8 root,
12 by 12 middle, and 16 by 16 finest plane-EB hierarchy. The split run writes
after its first accepted interval and must stop strictly before final time.
The checkpoint must carry the dedicated three-level magic and a complete end
marker before restart is attempted.

Restart must recover the accepted time and all three private field candidates,
advance to the same final time, and emit the expected 64, 144, and 256 rows.
Every root, middle, and finest CSV value must agree with the uninterrupted run
within a `3e-10` scale-aware tolerance. The focused reference,
checkpoint-stop, restart, and comparison gates must pass before all 187 tests
in GNU Fortran Debug and Release configurations.

## 0.104.0 dynamic three-level finest-patch gates

The public hotspot case begins with a 12 by 12 root, a fixed ratio-two
20 by 20 middle rectangle, and an active 8 by 8 finest seed. Initialization
must detect EB-aware temperature-gradient tags inside the qualified two-cell
middle margin and commit at least one topology change. Covered-side tags are
excluded, producing a deterministic 22 by 28 finest rectangle that crosses
the plane embedded boundary.

Root, middle, and finest outputs must reach `1e-7`, contain finite fields and
regular, cut, and covered cells, retain positive active thermodynamics, and
preserve species closure within `8e-12`. The focused unit, public run, and
output-structure gates must pass before all 189 tests in GNU Fortran Debug and
Release configurations.

## 0.105.0 dynamic three-level checkpoint/restart gates

The public split-run case begins from the same 12 by 12 root, fixed 20 by 20
middle rectangle, and configured 8 by 8 finest seed. Initialization must first
commit the EB-aware 22 by 28 finest topology. The checkpoint-stop run writes
after its first accepted interval, records that actual rectangle and at least
one committed regrid, and stops strictly before `2e-7`.

Restart must reconstruct the stored 22 by 28 finest field rather than the
configured seed, restore accepted step and regrid accounting, and continue
the same cadence to final time. The checkpoint checker requires the distinct
dynamic magic, exact regrid controls, stored non-seed topology, regrid count,
and terminal marker. Every root, middle, and finest CSV field must agree with
the uninterrupted run within a `3e-10` scale-aware tolerance. These four
focused gates must pass before all 193 tests in GNU Fortran Debug and Release
configurations.

## 0.106.0 single-level reactive EB transport gates

The low-level gate initializes a nonuniform temperature hotspot on a plane-EB
mesh, chooses a positive mixture transport stability limit, and advances
thermal conduction through the EB SSPRK2/StateRedist transaction. The active
temperature span must decrease while every fluid-volume-weighted conserved
integral remains fixed within an `8e-12` scale-aware tolerance. Covered cells
must remain bitwise unchanged, the species limiter must remain inactive for
the conduction-only case, and a negative interval must roll back exactly.

The public application runs matching inert and conducting 16 by 16 plane-EB
hotspots to `2e-7`. Both outputs must retain finite fields, regular/cut/covered
classes, positive active thermodynamics, and species closure within `8e-12`.
The conducting active-temperature span must be measurably smaller than the
inert result. A separate uniform-field gate exercises viscosity, conduction,
species diffusion, and barodiffusion together and requires state preservation
and composite conservation. These focused gates must pass before all 197 tests
in GNU Fortran Debug and Release configurations.

## 0.107.0 two-level reactive EB AMR transport gates

The low-level hierarchy gate initializes a nonuniform hotspot on nested
coarse/fine plane-EB meshes. It selects positive transport limits on both
levels, advances fine transport with ratio subcycling and time-interpolated
coarse exterior states, and requires the fine active-temperature span to
decrease. Diffusive reflux and reactive average-down must retain every
composite conserved integral within a `2e-10` scale-aware tolerance, preserve
covered cells bitwise, and roll back both levels for an invalid interval.

The driver gate must reject transport when its database is omitted, accept the
same request with the matching seven-species database, reach final time under
the combined coarse/fine parabolic limit, keep an active fine patch, and retain
composite conservation. The public application runs matching inert and
conducting 12 by 12 root plus 20 by 20 fine-patch hotspots to `2e-7`. All four
CSVs must contain finite positive active thermodynamics and species closure
within `8e-12`; the conducting fine-level temperature span must be measurably
smaller. These focused gates must pass before all 201 tests in GNU Fortran
Debug and Release configurations.

## 0.108.0 three-level reactive EB AMR transport gates

The low-level hierarchy gate initializes a nonuniform hotspot on nested root,
middle, and finest plane-EB meshes. It selects a positive transport limit on
every level, converts both child limits to the root interval, and advances
middle and finest transport with nested ratio subcycling and time-interpolated
parent exterior states. Independent diffusive registers must reflux each
interface in deepest-first order. Every three-level composite conserved
integral must remain within a `3e-10` scale-aware tolerance, covered cells must
remain bitwise unchanged, and an invalid interval must roll back all levels.

The driver gate must reject an omitted transport database and accept the same
request with the matching seven-species database. The public application runs
matching inert and conducting 8 by 8 root, 12 by 12 middle, and 16 by 16
finest hotspots to `2e-7`. All six CSVs must contain finite positive active
thermodynamics and species closure within `8e-12`; the conducting composite
hierarchy temperature span must be measurably smaller. These focused gates must pass
before all 205 tests in GNU Fortran Debug and Release configurations.

## 0.109.0 multipatch reactive EB AMR transport gates

The low-level gate initializes two separated ratio-two fine patches over a
plane-EB root and imposes a double temperature hotspot across the composite
hierarchy. The coarse transport stage must run once, while each child owns an
independent diffusive register and fine subcycle. A set-wide EB-cut closure
must retain every composite conserved integral within a `5e-10` scale-aware
tolerance. Covered cells remain bitwise unchanged, the conduction-only
species limiter remains inactive, the hierarchy temperature span decreases,
and an invalid interval rolls back the root and every child exactly.

The public application runs matching inert and conducting 14 by 14 root plus
two 10 by 10 sibling-patch double hotspots to `2e-7`. All six CSVs must retain
finite positive active thermodynamics, identical EB classification, and
species closure within `8e-12`. The conducting composite hierarchy
temperature span must be measurably smaller. These focused gates must pass
before all 208 tests in GNU Fortran Debug and Release configurations.

## 0.110.0 MPI reactive EB AMR ownership gates

The MPI gate constructs the qualified 14 by 14 plane-EB root and two separated
ratio-two 10 by 10 children independently on every rank. A collective
topology check must accept identical descriptors. The deterministic
parabolic-work schedule must account for 196 root cells plus two children at
`100 * r^2` work each, give every tested rank a root tile, and assign every
entity exactly once.

Each owner writes a distinct but physically valid scaled reactive payload
while nonowners retain stale finite replicas. Synchronization must recover the
exact owner state and temperature for every root tile and child on every rank.
An out-of-range owner must reject without changing outputs, a rank-dependent
work exponent must reject collectively, and exponent three must reject on all
ranks. The executable must pass with OpenMPI at one, two, four, and eight
ranks in GNU Fortran Release and bounds/FPE-checked Debug builds before the
complete 208-test serial regression in each MPI configuration.

## 0.111.0 owner-only MPI reactive EB AMR chemistry gates

The ownership executable builds a serial active-cell chemistry reference for
the same plane-EB root and two separated ratio-two children. Its MPI path must
call the reactor exactly once on the exclusive owner of every root tile and
child patch. Each rank's committed call count must equal its owner-map entity
count, and the global sum must equal the complete entity count.

After each accepted owner transaction, the owner state and EOS-recovered
temperature are broadcast. A fine-to-root average-down then must match the
serial patch-set reference to `5e-13` of the field scale for the root and both
children, while at least one reactive species changes measurably. Corrupting
the density on the owner of the last root tile must reject collectively after
earlier entities have advanced and must leave every rank's original root and
child fields bitwise unchanged. The gate runs at one, two, four, and eight
ranks in GNU Fortran Release and bounds/FPE-checked Debug configurations.

## 0.112.0 owner-only MPI reactive EB AMR hydro gates

The MPI path must reproduce the serial multipatch EB hydro transaction for a
nonmatching root plus two separated ratio-two children. The exclusive root
physics owner advances the complete weighted-StateRedist level exactly once.
Each child owner must execute both fine substeps, accumulate the matching
coarse and fine time-integrated fluxes, reflux the authoritative root
candidate in child order, and publish its corrected child. The root owner then
performs one patch-set average-down.

Each rank's committed advance count must equal one when it owns the root
physics interval plus the refinement ratio of every locally owned child. The
global count must be five for the qualified hierarchy. Root and child state
and temperature fields must match the serial transaction within `8e-12` of
their field scales and must change measurably from the nonuniform input.

The owner of the final child receives a finite but EOS-invalid density field.
The root and earlier child candidate work must complete before that fine
advance rejects collectively. Every rank must retain its original root and
child payloads bitwise and report zero committed advances. The gate runs at
one, two, four, and eight ranks in GNU Fortran Release and
bounds/FPE-checked Debug configurations.

## 0.113.0 owner-only MPI reactive EB AMR transport gates

The MPI path must reproduce the serial multipatch SSPRK2 molecular-transport
transaction for a plane-EB root and two separated ratio-two children with
different density, temperature, and velocity states. Before advancement,
nonowners hold deliberately stale finite replicas; owner synchronization must
still recover the serial start hierarchy.

In each Euler stage the exclusive root physics owner advances the complete EB
level once. Each child owner must execute both fine substeps, accumulate its
coarse and fine time-integrated diffusive fluxes, reflux in child order, and
publish the corrected child and root. The root owner performs the set-wide
average-down and the EB-cut composite conservation closure. Across two SSPRK2
stages, the qualified hierarchy must report exactly ten committed Euler-level
advances. State, temperature, and limiter minima must match the serial
transaction within their scale-aware tolerances and the nonuniform hierarchy
must change measurably.

The final child owner receives a finite but EOS-invalid density field after
the root and earlier child candidates can execute. The transaction must reject
communicator-wide, leave every rank's caller-owned root and child fields
bitwise unchanged, and report zero committed advances. The gate runs at one,
two, four, and eight ranks in GNU Fortran Release and bounds/FPE-checked Debug
configurations.

## 0.114.0 owner-only MPI reactive EB AMR full-physics gates

The MPI transaction must match the existing serial multipatch `R-T-H-T-R`
composition on the qualified plane-EB root and two ratio-two children. The
uniform reactive start field isolates composition and ownership: chemistry
provides a measurable state change while the dedicated 0.112 and 0.113 gates
retain the nonuniform hydro and transport coverage.

Across the successful outer interval, the communicator must report exactly
two chemistry calls for every root tile and child, five hydro level advances,
and twenty transport Euler-level advances. Root and child fields and the
transport limiter minimum must match the serial transaction within `5e-11`
field-scale and `3e-13` limiter tolerances.

A second transaction uses valid chemistry, transport, and state controls but
an invalid Riemann solver. The first chemistry and transport prefixes must
execute before hydro rejects. The caller's root and child fields must remain
bitwise unchanged and all three published counters must remain zero on every
rank. The gate runs at one, two, four, and eight ranks in GNU Fortran Release
and bounds/FPE-checked Debug configurations.

## 0.115.0 sparse MPI reactive EB AMR storage gates

Starting from deliberately distinct owner payloads and stale nonowner
replicas, scatter must allocate state and temperature only for locally owned
root tiles and fine children. Each rank's stored-value count must equal its
owner-map cell count times the reactive state-plus-temperature width, and the
communicator sum must equal one complete root plus both children exactly once.

Materialization must reproduce the established owner-authoritative root and
child synchronization bitwise at one, two, four, and eight ranks. Removing one
rank-zero-owned root payload must reject collectively before output commit;
every returned fallback root and child field must remain bitwise unchanged.
The gate runs in GNU Fortran Release and bounds/FPE-checked Debug
configurations before the complete 208-test regression.

## 0.116.0 direct sparse MPI reactive EB AMR chemistry gates

The sparse chemistry path must react only locally allocated root tiles and
children, with covered cells masked, and report exactly the distribution's
local entity count. The communicator sum must equal all root tiles and both
children once. After temporary synchronization and sparse re-scatter,
materialization must match the serial patch-set chemistry root and child fields
bitwise at one, two, four, and eight ranks.

For rollback, the last child owner receives a finite negative-density payload.
Root and earlier child reactors may execute, but collective rejection must
restore every locally allocated sparse state and temperature bitwise and report
zero committed calls. All rank-dependent comparisons must aggregate locally
before entering collectives. The gate runs in GNU Fortran Release and
bounds/FPE-checked Debug configurations before the complete 208-test
regression.

## 0.117.0 direct sparse MPI reactive EB AMR average-down gates

The direct sparse chemistry result must retain the exact local sparse value
count and materialize bitwise-identically to serial patch-set chemistry for
root state, root temperature, child state, and child temperature at one, two,
four, and eight ranks. This exercises child-owner volume-weighted restriction,
coarse-footprint broadcast, covered-cell preservation, and root-owner EOS
temperature recovery without a complete temporary hierarchy.

The rejection gate gives the final child owner finite negative-density state.
The independent sparse average-down call must reject collectively after any
earlier child restrictions, preserve every local root and child allocation
bitwise, and retain its exact local value count. The gate runs in GNU Fortran
Release and bounds/FPE-checked Debug configurations before the complete
208-test regression.

## 0.118.0 sparse MPI reactive EB AMR full-physics gates

Starting from an owner-only sparse hierarchy, the complete `R-T-H-T-R`
transaction must retain the exact local sparse value count and match the
serial multipatch full-physics root and child state and temperature within
`5e-11` field scale. The transport limiter minimum must match within `3e-13`.
At one, two, four, and eight ranks, local and communicator counts must equal
two chemistry calls per root tile and child, one root plus ratio-subcycled
child hydro update, and four root plus ratio-subcycled child transport Euler
updates.

The rollback transaction uses valid chemistry and transport controls but an
invalid hydro solver. Chemistry and the first transport stage may execute in
the candidate compatibility window, but the sparse caller state must remain
bitwise unchanged, all three published operator counts must remain zero, and
the limiter fallback must remain exactly one. The gate runs in GNU Fortran
Release and bounds/FPE-checked Debug configurations before the complete
208-test regression.

## 0.119.0 direct sparse MPI reactive EB AMR hydro gates

Starting from sparse root tiles and owner-only fine children, direct hydro must
retain each rank's exact stored-value count. Materialization only after commit
must match serial multipatch hydro root and child state and temperature within
`8e-12` field scale at one, two, four, and eight ranks. Local and communicator
counts must equal one root-owner update plus each child refinement ratio on its
exclusive owner.

The rejection gate places finite negative density on the final child owner.
The root and any earlier children may advance first, but every local sparse
root and child field must remain bitwise unchanged and the published advance
count must remain zero. The gate runs in GNU Fortran Release and
bounds/FPE-checked Debug configurations before the complete 208-test
regression.

## 0.120.0 direct sparse MPI reactive EB AMR transport gates

Starting from sparse root tiles and owner-only fine children, direct SSPRK2
transport must retain each rank's exact stored-value count. Materialization
only after commit must match serial multipatch transport root and child state
and temperature within `2e-11` field scale and match the limiter minimum within
`2e-13` at one, two, four, and eight ranks. Local and communicator counts must
equal two root-owner Euler updates plus two times each child refinement ratio
on its exclusive owner.

The rejection gate places finite negative density on the final child owner.
The first root stage and any earlier child candidates may execute, but every
local sparse root and child field must remain bitwise unchanged, the published
Euler count must remain zero, and the limiter fallback must remain exactly one.
The gate runs in GNU Fortran Release and bounds/FPE-checked Debug configurations
before the complete 208-test regression.

## 0.121.0 end-to-end sparse MPI reactive EB AMR full-physics gates

Starting from owner-only sparse root tiles and fine children, the complete
`R-T-H-T-R` transaction must retain each rank's exact stored-value count and
match serial multipatch full physics for root and child state and temperature
within `5e-11` field scale. The transport limiter minimum must match within
`3e-13`. Local and communicator counts must match two sparse chemistry
half-steps, one sparse hydro interval, and two sparse SSPRK2 transport
half-steps on the exclusive physics owners at one, two, four, and eight ranks.

The rejection gate uses valid chemistry and transport controls with an invalid
hydro solver. The first sparse chemistry and transport candidates may succeed,
but the caller's local sparse fields must remain bitwise unchanged, all three
published operator counts must remain zero, and the limiter fallback must
remain exactly one. The gate runs in GNU Fortran Release and bounds/FPE-checked
Debug configurations before the complete 208-test regression.

## 0.122.0 targeted sparse MPI reactive EB AMR average-down gates

For every child, derive the unique root tile owners whose row bands intersect
its coarse footprint, excluding the child owner itself. The child owner must
report exactly one point-to-point restriction transfer per remaining owner;
the communicator sum must equal the independently computed recipient count at
one, two, four, and eight ranks. With one rank the count must be zero.

The chemistry path that consumes this average-down must retain its exact sparse
stored-value count and reproduce serial root and child state and temperature
bitwise. A direct average-down rejection from finite negative child density may
send candidate restrictions, but every local sparse field must remain bitwise
unchanged and the published transfer count must remain zero. The gate runs in
GNU Fortran Release and bounds/FPE-checked Debug configurations before the
complete 208-test regression.

## 0.123.0 targeted direct sparse MPI reactive EB AMR hydro gates

For a root split into one row tile per available rank, the expected hydro
payload count is two transfers per non-root-owner tile, one root bundle per
distinct remote child owner, and two correction transfers per remote child.
Each rank's reported sends and the communicator sum must match that independent
formula at one, two, four, and eight ranks; the one-rank count must be zero.

The direct sparse hydro result must retain every rank's exact stored-value
count and match serial multipatch hydro root and child state and temperature
within `8e-12` field scale. Finite negative density on the final child may
follow successful root and earlier-child candidates, but all local sparse
fields must remain bitwise unchanged and both published advance and transfer
counts must remain zero. The gate runs in GNU Fortran Release and
bounds/FPE-checked Debug configurations before the complete 208-test
regression.

## 0.124.0 targeted direct sparse MPI reactive EB AMR transport gates

For each of the two Euler stages, the expected payload count is one gather and
one final scatter per non-root-owner tile, one bundle per distinct remote child
owner, and two correction transfers per remote child. The final SSPRK2 blend
adds two gathers and one scatter per non-root-owner tile. Each rank's reported
sends and the communicator sum must match this independent formula at one,
two, four, and eight ranks; the one-rank count must be zero.

The direct sparse SSPRK2 result must retain every rank's exact stored-value
count and match serial multipatch root and child state and temperature within
`2e-11` field scale. Its limiter minimum must match within `2e-13`. Finite
negative density on the final child may follow successful root and
earlier-child candidates, but all local sparse fields must remain bitwise
unchanged, the limiter fallback must remain one, and both published Euler and
transfer counts must remain zero. The gate runs in GNU Fortran Release and
bounds/FPE-checked Debug configurations before the complete 208-test
regression.

## 0.125.0 sparse owner-local MPI reactive EB AMR timestep gates

Build a serial reference from the root and every fine child's EB hyperbolic
CFL limit plus the enabled molecular-transport stability limit. Scale every
fine-level limit by its refinement ratio before taking the hierarchy minimum.
The sparse owner-only result must match that reference within 64 machine
epsilons at one, two, four, and eight ranks.

The expected communication count is one root gather per non-root-owner tile.
Each rank's reported sends and the communicator sum must match this independent
formula, with a zero count at one rank. A finite negative density on the final
child must reject collectively after any root gather, preserve all sparse
fields bitwise, and publish zero dt and zero transfers. The gate runs in GNU
Fortran Release and bounds/FPE-checked Debug configurations before the complete
208-test regression.

## 0.126.0 public sparse MPI reactive EB AMR time-loop gates

Choose a target time equal to `1.25` times the initial serial full-physics
stable interval so at least two dynamically selected steps are required and
the final one is clipped. An independent serial patch-set loop recomputes the
combined EB hydro/transport limit after each accepted state. The sparse public
loop must reach the target exactly and match the serial step count, minimum dt,
limiter minimum, root state and temperature, and every child state and
temperature at one, two, four, and eight ranks.

Local and communicator-summed chemistry, hydro, and transport counts must equal
the independently expected per-step owner counts multiplied by the committed
step count. Timestep root-gather traffic must obey the same multiplication. A
second run with a total-step limit of one must reject before the next interval,
retain exactly the first committed state, clock, minimum dt, and diagnostics,
and leave the sparse hierarchy valid. The gate runs in GNU Fortran Release and
bounds/FPE-checked Debug configurations before the complete 208-test
regression.

## 0.127.0 transactional sparse MPI reactive EB AMR regrid gates

Start from two separated ratio-two child patches and give each child a distinct
but thermodynamically consistent conserved-state scale. Build an independent
serial reference that replaces them with one shifted and resized child. The
sparse regrid must reproduce the serial averaged root, retained fine overlap,
newly prolonged fine cells, temperatures, and ordered topology exactly at one,
two, four, and eight ranks.

After ownership is recomputed, each rank's stored values must equal `nvar+1`
times its assigned root/child cell count, and the communicator sum must contain
exactly one numerical copy of the rebuilt hierarchy. A refinement ratio of one
must reject before materialization and preserve the old distribution, child
owners, sparse fields, template, and `changed=false` result. The gate runs in
GNU Fortran Release and bounds/FPE-checked Debug configurations before the
complete 208-test regression.

## 0.128.0 scheduled tag-driven sparse MPI reactive EB AMR regrid gates

Start from two separated ratio-two children and a thermodynamically consistent
root temperature hotspot. Advance an independent serial full-physics loop to
`1.01` times the initial stable interval, evaluating temperature tags after
accepted step two and applying the serial multipatch regrid. The public sparse
loop must match its exact step, evaluation, and topology-change counts; minimum
dt and limiter; ordered topology; and root/child state and temperature within
`2e-10` field scale at one, two, four, and eight ranks.

The single scheduled plan must send one root payload per non-root-owner tile,
while timestep selection sends the same independently counted payloads once per
accepted step. A second run uses a geometry builder that rejects the first
tagged plan. The otherwise valid physics candidate, time, step, distribution,
sparse fields, template, regrid counters, and published timestep/regrid traffic
must all remain at their pre-step values. The gate runs in GNU Fortran Release
and bounds/FPE-checked Debug configurations before the complete 208-test
regression.

## 0.129.0 direct sparse MPI reactive EB AMR regrid gates

Reuse both the explicit shifted/resized-child case and the scheduled hotspot
case. Their materialized root, retained fine overlap, newly PCM-prolonged cells,
recovered temperatures, ordered topology, and one-copy storage must continue to
match the independent serial references at one, two, four, and eight ranks.

For each old child, restriction sends must equal its distinct intersecting root
owners excluding its own owner. For each distinct new child owner, PCM root
assembly sends must equal root tiles owned elsewhere. Each nonempty same-ratio
old/new overlap rectangle whose owner changes must produce exactly one send.
Local sender counts and communicator sums are checked independently for the
explicit API and through the scheduled tagged public clock. Invalid controls
and a valid but overlap-inconsistent fine geometry must publish zero
direct-regrid traffic and preserve all caller state, including when rejection
occurs after restriction and PCM staging. A rejected scheduled geometry
callback obeys the same publication rule. Static inspection must show that the
explicit transaction no longer calls the all-rank materialization helper. The
gate runs in GNU Fortran Release and bounds/FPE-checked Debug configurations
before the complete 208-test regression.

## 0.130.0 root-only sparse MPI reactive EB AMR materialization gates

Scatter a valid two-child reactive EB patch set into exclusive sparse root
tiles and child payloads. Gather it to the deterministic root writer and require
bitwise equality with the established owner-authoritative materialization for
the complete root state, root temperature, and every child field. Only the
writer rank may allocate those complete outputs; all non-writer arrays remain
unallocated and their returned patch sets remain empty.

The expected traffic is one send for every root tile and child whose owner is
not the selected writer. Each sender-local count and the communicator sum must
match that independent owner-map formula at one, two, four, and eight ranks.
Removing one local sparse root field must reject collectively, deallocate all
outputs, return an empty patch set, and publish zero traffic. The gate runs in
GNU Fortran Release and bounds/FPE-checked Debug configurations before the
complete 208-test regression.

## 0.131.0 sparse MPI reactive EB AMR I/O gates

Write the exclusively owned two-child sparse hierarchy through the root-only
checkpoint adapter. Read the resulting file with the established serial
multipatch reader and require exact conserved-field parity, EOS-consistent
temperature parity, matching topology, clock, step, regrid, minimum-timestep,
and base-density metadata. The successful sender-local and communicator-summed
counts must equal the independent remote-entity formula at one, two, four, and
eight ranks.

Write the same hierarchy through the root/child CSV adapter and require one
nonempty root file plus one nonempty deterministically named file per child.
Checkpoint and CSV paths under a nonexistent directory must fail collectively
and publish zero traffic. Both adapters run in GNU Fortran Release and
bounds/FPE-checked Debug configurations before the complete 208-test serial
regression.

## 0.132.0 root-only sparse MPI reactive EB AMR restart gates

First gather the exclusively owned hierarchy to one root, then scatter those
root-only arrays directly back to the configured owners. Every owner-local root
tile and child field must be bitwise identical to the original sparse payload,
and the stored-value count must remain exact. The selected root reports one send
per entity owned elsewhere; all other ranks report zero, and the communicator
sum must match that independent owner-map formula.

Next write and read the formatted checkpoint through the public sparse I/O
adapters. Require the same owner-local field parity, transfer formula, and exact
clock metadata at one, two, four, and eight ranks. A missing checkpoint must
return an empty sparse set, zero clock metadata, and zero traffic collectively.
The gate runs in GNU Fortran Release and bounds/FPE-checked Debug configurations
before the complete 208-test serial regression.

## 0.133.0 geometry-only sparse MPI reactive EB AMR restart gates

Extract a descriptor containing only each child's EB geometry and patch box,
then perform the direct root-to-owner scatter without passing replicated child
state or temperature fields. Every owner-local root tile and child must remain
bitwise identical to the original sparse payload, with the same exact local and
communicator-summed remote-entity transfer formula at one, two, four, and eight
ranks.

Read the formatted checkpoint on the selected root through the geometry-only
API and require the same field, stored-value, and clock-metadata parity. A
descriptor with invalid EB volume fractions must reject collectively before
point-to-point traffic and publish an empty sparse set. A missing checkpoint
must preserve the established empty-state, zero-metadata, zero-traffic result.
The gate runs in GNU Fortran Release and bounds/FPE-checked Debug configurations
before the complete 208-test serial regression.

## 0.134.0 arbitrary-depth geometry-only reactive EB topology gates

Build two separated root children, refine both branches once, and refine one
of those grandchildren again. The resulting four-level topology must validate
its ordered parent links, flattened child offsets and indices, per-level patch
counts, refinement ratios, full EB geometry consistency, and sibling
separation.

Rebuild the accepted three-level tree to the four-level plan and require one
committed topology change. Repeating the identical plan must report a no-op.
Changing the deepest child to a nonexistent parent must reject without changing
the accepted four-level topology. The gate runs inside the established EB
multilevel unit in GNU Fortran Release and bounds/FPE-checked Debug before the
complete 208-test serial regression.

## 0.135.0 arbitrary-depth reactive EB state-migration gates

Initialize reactive fields on the established three-level, two-branch tree,
give both deeper branches distinct valid conserved states, and require the
composite conserved vector to remain unchanged after deepest-first
synchronization. Rebuild to the four-level plan and require parent-first PCM
initialization, same-resolution overlap retention, recovered positive finite
temperatures, and conservation of every state component.

Move the deepest rectangle by one parent cell after giving its fine cells a
nonuniform but EOS-valid state. The physical intersection must retain the old
fine values under the expected index shift, newly exposed cells must come from
the updated parent, and the complete composite integral must remain unchanged.

Repeating the four-level plan must be an exact numerical no-op. A deepest child
that names a nonexistent parent must reject with the accepted topology, state,
and temperature bitwise unchanged. The gate runs inside the established EB
multilevel unit in GNU Fortran Release and bounds/FPE-checked Debug before the
complete 208-test serial regression.

## 0.136.0 MPI owner-tiled reactive EB root-hydro gates

Replace the replicated owner path's selected-rank root advance with one
finite-halo band advance per root tile. At one, two, four, and eight ranks,
require the assembled root state, recovered temperature, all Cartesian/EB face
flux effects, every fine child, reflux, and final average-down to match the
established serial multipatch transaction within the qualified tolerances.

Count exactly one root advance on each tile owner rather than one communicator
advance on the first tile owner. Independently sum the actual halo-band cell
counts and require the exact distribution-derived value; above one rank it must
be smaller than computing the complete root independently on every rank. A late
child rejection must publish zero advances and zero root work. Repeat the same
owner accounting and serial parity through the complete chemistry-transport-
hydro split. Run the gates with OpenMPI at one, two, four, and eight ranks in
Release and bounds/FPE-checked Debug configurations before the 208-test serial
regression.

## 0.137.0 sparse MPI owner-tiled reactive EB root-hydro gates

Starting from exclusively owned sparse root tiles, construct each target
tile's six-row EB band from local rows plus point-to-point fragments sent only
by intersecting source owners. Require exactly one bounded-band advance per
root tile owner and the exact distribution-derived computed-cell count. Above
one rank, total sparse root work must be smaller than independently advancing
the full root on every rank.

Route each tile's owned input, updated state, temperature, x-flux rows, and
uniquely owned y-faces to the root owner, then retain the established targeted
child bundle, correction round trips, and final row scatter. Count every
remote halo and result payload exactly. At one, two, four, and eight ranks,
materialized root and child fields must match the serial multipatch hydro
transaction within `8e-12` field scale. Repeat tile-owner call accounting
through sparse `R-T-H-T-R` and the public time loop. A late child rejection
must leave all sparse fields bitwise unchanged and publish zero advances,
computed cells, and transfers. Run the gates in OpenMPI Release and
bounds/FPE-checked Debug before the complete serial regression.

## 0.138.0 sparse MPI owner-local EB timestep gates

Evaluate the hyperbolic CFL and optional molecular-transport stability limit
directly on every exclusively owned root tile and fine child. Root tile calls
must use the exact extracted EB row-band geometry; fine limits must retain
their refinement-ratio conversion to a coarse interval. A communicator minimum
must match the serial complete-hierarchy timestep within `64 epsilon` at one,
two, four, and eight ranks while reporting exactly zero root transfers.

The public clipped time loop and scheduled-regrid loop must also publish zero
timestep root traffic on every accepted step. Regrid-planning traffic remains
counted independently. A finite invalid child state must reject collectively
with zero dt and zero transfers while preserving every sparse field bitwise.
Run the gates in OpenMPI Release and bounds/FPE-checked Debug before the
complete serial regression.

## 0.139.0 sparse MPI owner-local SSPRK2 root-blend gates

After both sparse transport Euler stages succeed, blend every exclusively
owned root tile directly from its local interval-start and second-Euler
candidates. Recover temperature with the exact extracted EB row-band geometry.
At one, two, four, and eight ranks, require root and child state and temperature
to retain the established serial transport tolerances, and require the limiter
minimum to retain its established parity tolerance.

For every remote root tile, require exactly four root transfers across the two
Euler stages: one gather and one scatter per stage. The final blend must add no
start gather, second-Euler gather, or blended-state scatter. Keep child bundles,
reflux correction traffic, distributed cut-interface closure, stored-value
counts, and late-failure rollback unchanged. Run in OpenMPI Release and
bounds/FPE-checked Debug before the complete serial regression.

## 0.140.0 sparse MPI owner-tiled SSPRK2 root-Euler gates

For each Euler stage and root tile, assemble the exact six-row transport guard
from local rows plus direct fragments from every intersecting remote source
owner. Execute one EB transport-flux/StateRedist band advance on the target tile
owner, then route its owned start, result, temperature, x-flux, and unique
y-face rows to the root physics owner. Require exact distribution-derived halo,
result, child-bundle, correction, and scatter transfer counts.

Across a complete SSPRK2 call, require two tile advances per root tile and the
exact sum of both stages' band-cell counts. On the small qualification mesh,
the periodic-edge safeguard may cover the full root through two ranks; above
two ranks total work must remain below independently advancing the full root on
every rank. Materialized root and child fields must retain the established
`2e-11` transport tolerance and the limiter minimum its `2e-13` tolerance at
one, two, four, and eight ranks. Repeat owner accounting through sparse
`R-T-H-T-R` and the public time loop. A late failure must preserve every sparse
field bitwise and publish zero work and traffic. Run Release and bounds/FPE-
checked Debug before the complete serial regression.

## 0.141.0 sparse MPI periodic-edge cyclic-band gates

Build a dedicated 14-by-21, root-only periodic-y EB transport case. At four
and eight ranks, require each boundary target to use two increasing global
source-row fragments whose combined band is smaller than the complete root.
The band must retain the six-row transport/StateRedist dependency footprint
plus one seam-isolation row, while one and two ranks exercise the complete-root
fallback.

Require exact owner advances, point-to-point fragment and result transfers,
and computed-band-cell totals derived independently in the test. Materialize
the sparse result and compare state and temperature with a serial SSPRK2
transport reference at the established `2e-11` scaled tolerance; compare the
redistribution limiter minimum at `2e-13`. Repeat the existing one-, two-,
four-, and eight-rank full-physics, public time-loop, and scheduled-regrid
gates in OpenMPI Release and bounds/FPE-checked Debug, followed by the complete
serial regression.

## 0.142.0 compact EB flux-register and MPI correction gates

Initialize an interior EB flux register and require its correction array to
retain absolute coarse indices while storing only the fine-patch rectangle
expanded by one cell. Re-run aperture-matched coarse/fine cancellation,
cut-cell extensive conservation, fine-covered redistribution, reactive
temperature recovery, reset, nonfinite accumulation rollback, and nonphysical
reflux rollback without relaxing their established tolerances.

For sparse MPI transport, replace each remote child's full-root correction
round trip with the patch rectangle expanded by two coarse cells. The second
cell covers every cardinal or diagonal recipient of a cut-cell correction.
Require unchanged exact message counts, owner work counts, limiter minima, and
serial root/child field tolerances at one, two, four, and eight ranks. Repeat
full physics, public time-loop, scheduled-regrid, and late-failure rollback in
OpenMPI Release and bounds/FPE-checked Debug before the complete serial
regression.

## 0.143.0 compact sparse child transport-context gates

Extract the four-edge start/end coarse context and require its reconstructed
exterior state and temperature to match the complete-root builder exactly.
For every qualification child, require the context-plus-coarse-register value
count to be positive and smaller than the former complete root
start/end/temperature/x-y-flux bundle.

For a remote child and each transport Euler stage, require exactly three
messages: compact context and coarse mismatch from root to child, evolved fine
state plus accumulated mismatch from child to root, and corrected fine state
from root to child. Reflux remains ordered on the root owner. Retain exact
owner work accounting, limiter minima, serial root/child tolerances, public
time-loop accounting, and late-failure rollback at one, two, four, and eight
ranks in Release and bounds/FPE-checked Debug.

## 0.144.0 compact child-local reactive reflux gates

Run reactive reflux once through complete-root compatibility arrays and once
through a strictly smaller array passed to the global-index support entrypoint,
and require bitwise-identical coarse state and temperature on the supplied
support plus bitwise-identical complete fine state and temperature. Both
successful registers must reset; established nonphysical rollback remains
unchanged.

For sparse MPI transport, require the combined exterior/register/patch-plus-
two context payload to be smaller than the former complete root bundle. A
remote child must use exactly two messages per Euler stage: context/support
from root to child and corrected support from child to root. Retain exact
root/child field tolerances, limiter minima, owner work, public-clock traffic,
scheduled regridding, and late-failure rollback at one, two, four, and eight
ranks in Release and bounds/FPE-checked Debug.

## 0.145.0 compact coarse interface-flux gates

Accumulate one coarse flux register from complete root x/y flux arrays and a
second from globally indexed rectangles containing only the four coarse/fine
interface face ranges. Require bitwise-identical corrections and require the
compact value count to be strictly smaller. Omit one active interface face and
require rejection with a bitwise-unchanged register.

Run sparse MPI transport through the compact entrypoint while retaining exact
root/child fields, register reset, limiter minima, owner work, traffic, public
clock, scheduled-regrid, and rollback gates at one, two, four, and eight ranks
in Release and bounds/FPE-checked Debug.

## 0.146.0 direct root-tile coarse-flux routing gates

Retain each root tile's x-flux rows and uniquely owned y-faces, then assemble
each child's compact interface rectangles from only intersecting tile owners.
Require complete receiver coverage, finite values, and a combined direct-flux
plus state-context value count below the former complete root bundle.

Remove the coarse register from the root-to-child state context. Derive the
exact point-to-point count from every remote tile/child intersection plus the
existing context and corrected-support messages. Require unchanged serial
root/child fields, register reset, limiter minima, owner work, public clock,
scheduled regridding, and late-failure rollback at one, two, four, and eight
ranks in Release and bounds/FPE-checked Debug.

## 0.147.0 compact exterior state-context gates

Extract one child exterior context from complete root start/end state and
temperature arrays and a second from a strictly smaller globally indexed
patch-plus-one support. Reconstruct at an interior coarse-time fraction and
require bitwise-identical state and temperature on all four fine edges.

Require incomplete, out-of-root, and nonfinite support to reject without a
valid published context. Retain all established serial reactive EB AMR and
sparse MPI transport gates before qualifying the support API as the next
direct state-routing boundary.

## 0.148.0 direct root-tile state/support routing gates

Retain stage-start, uncorrected stage-end, and current corrected root-tile
state/temperature. Assemble each child's patch-plus-two support directly from
intersecting owners and extract its exterior context on the child. Require the
resulting state plus interface-flux value count to remain below the former
complete root bundle.

Route the child-local reflux result back to every intersecting tile owner before
the next child begins, then commit final root tiles without root-owner scatter.
Derive exact point-to-point counts for state, flux, and correction fragments
across both SSPRK2 Euler stages. Retain serial root/child fields, limiter minima,
owner work, public clock, scheduled regridding, cut-boundary conservation, and
late-failure rollback at one, two, four, and eight ranks in Release and
bounds/FPE-checked Debug.

## 0.149.0 owner-local root transport result gates

Remove the transport Euler stage-start, stage-end, temperature, and x/y-flux
result message from every remote root tile to the root physics owner. Require
the exact point-to-point count to contain only finite-band halo exchanges and
direct state, flux, and correction fragments for intersecting children across
both SSPRK2 Euler stages. The root-only cyclic case must contain halo traffic
only.

For a cut interface, compute left/right physical-boundary contributions on
every local tile and lower/upper contributions only on the edge tiles, then
combine one finite `nvar` vector collectively. Retain serial root/child fields,
cut-boundary conservation, limiter minima, owner work, public clock, scheduled
regridding, and late-failure rollback at one, two, four, and eight ranks in
Release and bounds/FPE-checked Debug before the complete serial regression.

## 0.150.0 compact sparse hydro child-context gates

For each hydro child, extract the complete-root four-edge start/end context,
current patch-plus-two corrected state/temperature, and intersecting coarse
x/y flux rectangle. Require the combined payload to be positive and strictly
smaller than the former complete root start/end/temperature/flux bundle.

For every remote child require one packed context/support/flux message from the
root physics owner and one corrected-support message back. Remove the distinct-
child-owner complete bundle allowance from exact traffic. Retain serial
root/child field parity, deterministic overlapping-child correction, owner
work, public clock, scheduled regridding, conservation, and late-failure
rollback at one, two, four, and eight ranks in Release and bounds/FPE-checked
Debug before the complete serial regression.

## 0.151.0 direct hydro coarse-flux routing gates

Retain each hydro root tile's x-flux rows and unique y-faces, then assemble each
child's interface rectangles from only intersecting tile owners. Require full
receiver coverage, finite values, and a positive combined state-context plus
interface-flux value count smaller than the former complete root bundle.

Remove flux values from both the tile-to-root hydro result and root-to-child
state context. Derive exact point-to-point traffic from finite-band halos,
remote tile state results, remote final row scatters, remote child state and
correction messages, and one flux message per remote tile/child intersection.
Retain serial root/child field parity, deterministic overlapping-child
correction, owner work, public clock, scheduled regridding, conservation, and
late-failure rollback at one, two, four, and eight ranks in Release and bounds/
FPE-checked Debug before the complete serial regression.

## 0.152.0 owner-local hydro result gates

Retain stage-start, uncorrected stage-end, and current corrected hydro state and
temperature on root tile owners. Assemble each child's patch-plus-two state
support from intersecting tile owners, extract the exterior context on the
child, and return reflux corrections directly before assembling the next child.
Require the state-support plus interface-flux value count to remain positive
and smaller than the former complete root bundle.

Remove tile-to-root state results, complete root hydro result allocation,
root-owner correction merge, and final row scatter. Derive exact point-to-point
traffic from finite-band halos and direct state, flux, and correction fragments
only. Retain serial root/child fields, deterministic overlapping-child
correction, owner work, public clock, scheduled regridding, conservation, and
late-failure rollback at one, two, four, and eight ranks in Release and bounds/
FPE-checked Debug before the complete serial regression.

## 0.153.0 arbitrary-depth reactive EB patch-tree timestep gates

Evaluate a four-level, two-branch numerical EB tree with the established
single-node active-cell CFL kernel. Independently reduce every node-local limit
after multiplication by the cumulative refinement product and require the new
tree selector to match, skipping the fully covered branch as a nonconstraint.
Make the deepest node uniquely limiting and require the published root interval
to equal its local interval times all three refinement ratios.

Require the selector to leave every state and temperature value unchanged. A
nonfinite CFL must reject with zero output and the same read-only contract.
Retain all serial tests in GNU Fortran Release and bounds/FPE-checked Debug, then
retain the established OpenMPI one-, two-, four-, and eight-rank suite even
though this milestone adds no MPI communication.

## 0.154.0 arbitrary-depth reactive EB patch-tree hydro gates

Advance the existing four-level, two-branch EB numerical tree with PCM and the
qualified level kernel. Require committed per-level node counts `[1, 4, 8, 8]`,
positive finite temperatures, a changed solution, and conservation of density,
total energy, and every species in the stationary-boundary composite integral.

Represent the existing cut-interface fixed three-level hierarchy as a runtime
patch-tree chain. Require `[1, 2, 4]` scheduling and field/temperature agreement
with the qualified fixed-depth advance within the documented floating-point
tolerance. An invalid solver must reject with zero published counts and exact
tree rollback. Retain all 208 serial tests in GNU Fortran Release and bounds/
FPE-checked Debug, plus the established OpenMPI one-, two-, four-, and eight-
rank suite even though the new tree remains serial.

## 0.155.0 arbitrary-depth reactive EB patch-tree chemistry gates

Advance chemistry on the four-level, two-branch numerical tree and require one
active-mask chemistry call per runtime patch, giving per-level counts
`[1, 2, 2, 1]`. Compose `R-H-R` Strang splitting and require committed
chemistry counts `[2, 4, 4, 2]` plus hydro counts `[1, 4, 8, 8]`. Preserve
composite density and total energy, density/species closure, positive finite
thermodynamics, and measurable species activity.

Represent the existing fixed three-level chemistry/hydro hierarchy as a
runtime chain. Require chemistry counts `[2, 2, 2]`, hydro counts `[1, 2, 4]`,
and field/temperature agreement with the qualified fixed-depth Strang path. An
invalid solver after the first chemistry half-step must preserve every accepted
node exactly and publish zero chemistry and hydro counts. Retain all 208 serial
tests in GNU Fortran Release and bounds/FPE-checked Debug, plus the established
OpenMPI one-, two-, four-, and eight-rank suite even though the new tree remains
serial.

## 0.156.0 arbitrary-depth reactive EB patch-tree transport gates

Represent the qualified fixed three-level transport hierarchy as a runtime
tree and require SSPRK2 node counts `[2, 4, 8]`, positive limiter theta, and
state/temperature agreement with the fixed-depth implementation. Retain its
composite conservation, conduction smoothing, and covered-cell invariance
gates.

Build a separate four-level tree with one middle patch, two separated children,
and one deeper child. Require the actual branching SSPRK2 schedule
`[2, 4, 16, 16]`, changed state, positive limiter theta, valid thermodynamics,
and complete composite-vector conservation. A negative interval must preserve
every node exactly, return theta one, and publish zero counts. Retain all 208
serial tests in GNU Fortran Release and bounds/FPE-checked Debug, plus the
established OpenMPI one-, two-, four-, and eight-rank suite even though the new
tree remains serial.

## 0.157.0 arbitrary-depth reactive EB patch-tree full-physics gates

Compose chemistry, SSPRK2 transport, recursive hydro, SSPRK2 transport, and
chemistry on one private runtime tree. For the fixed three-level chain, require
chemistry counts `[2, 2, 2]`, transport Euler counts `[4, 8, 16]`, hydro counts
`[1, 2, 4]`, positive limiter theta, and field/temperature agreement with the
qualified fixed-depth `R-T-H-T-R` implementation.

On the four-level branching tree, require chemistry counts `[2, 2, 4, 2]`,
transport counts `[4, 8, 32, 32]`, hydro counts `[1, 2, 8, 8]`, composite
conservation, and valid thermodynamics. An invalid solver after the first
reaction and transport prefixes must restore every node exactly, publish zero
counts, and return theta one. Retain all 208 serial tests in GNU Fortran Release
and bounds/FPE-checked Debug, plus the established OpenMPI one-, two-, four-,
and eight-rank suite even though the new tree remains serial.

## 0.158.0 public arbitrary-depth reactive EB patch-tree time-loop gates

On the fixed three-level runtime chain, require the combined selector to equal
the independently computed minimum of the all-node hyperbolic limits and the
root-time-scaled transport limits. Advance to a target requiring two clipped
steps and require exact tree parity with an independently repeated sequence of
stable-step selection plus full-physics transactions.

Require exact final time, total and newly advanced step counts, minimum
accepted interval, minimum limiter theta, and accumulated per-level chemistry,
transport, and hydro schedules. Limit the same run to one step and require a
failed return that retains exactly that committed prefix. Reject an invalid
solver on the first step with exact tree/time/count rollback and neutral
timestep/theta outputs. Retain all 208 serial tests in GNU Fortran Release and
bounds/FPE-checked Debug, plus the established OpenMPI one-, two-, four-, and
eight-rank suite even though the new tree remains serial.

## 0.159.0 MPI arbitrary-depth EB patch-tree ownership gates

Build a four-level topology with patch counts `[1, 1, 2, 1]`. At one, two,
four, and eight ranks, require deterministic ownership, exact global totals of
five nodes and 352 allocated cells, and subcycle-squared weighted work 6016.
Require every owner to publish a distinct node value and every rank to receive
the exact complete tree, with global publication accounting equal to five.

Make one rank's candidate invalid and require collective rejection, zero
publication accounting, and exact rank-local rollback. Also provide unequal
subcycle exponents across ranks and require collective initialization
rejection; use an out-of-range exponent for the one-rank case. Retain the full
existing MPI gates and all 208 serial tests in GNU Fortran Release and bounds/
FPE-checked Debug.

## 0.160.0 MPI sparse arbitrary-depth EB patch-tree storage gates

Convert the accepted four-level branching tree into owner-only storage at one,
two, four, and eight ranks. Require each rank to allocate exactly its mapped
node and cell totals, with no fields on nonowners. Explicitly materialize the
sparse tree and require exact state and temperature parity with the accepted
replicated tree plus five global owner publications.

Rotate every owner by one rank when multiple ranks are present, migrate each
changed node directly, and require the global sender transfer count to equal
the number of ownership changes. Materialize again and require exact field
parity. Corrupt the new owner map on only rank zero and require collective
rejection, zero transfers, and value-for-value sparse rollback. Retain the
full established MPI gates and all 208 serial tests in GNU Fortran Release and
bounds/FPE-checked Debug.

## 0.161.0 MPI owner-local arbitrary-depth EB timestep gates

Reinitialize the accepted physical four-level branching tree into sparse
storage after rotating ownership. Enable hydro, viscosity, thermal conduction,
and species diffusion. Require each rank to evaluate only its owned active
nodes, the global evaluation sum to equal five, and the communicator minimum
to equal the complete serial patch-tree timestep exactly.

On multiple ranks, supply a different but individually valid hydro CFL on rank
zero and require collective rejection before node evaluation. On one rank, use
an invalid negative CFL for the same rejection path. In both cases require zero
timestep and zero local evaluation accounting. Retain the full established MPI
gates and all 208 serial tests in GNU Fortran Release and bounds/FPE-checked
Debug.

## 0.162.0 MPI owner-local arbitrary-depth EB chemistry gates

Advance the physical four-level branching tree for one chemistry interval
after rotating ownership. Require the global per-level owner advances to equal
`[1, 1, 2, 1]`, every parent/child owner difference to produce exactly one
direct child-state transfer, and shared-owner restriction to produce none.
Materialize only after the operation and require exact state and temperature
parity with the complete serial chemistry transaction.

On multiple ranks, supply a different but valid interval on rank zero and
require collective preflight rejection. On one rank, use a negative interval.
In both cases require exact sparse rollback, zero level advances, and zero
restriction transfers. Retain the full established MPI gates and all 208
serial tests in GNU Fortran Release and bounds/FPE-checked Debug.

## 0.163.0 MPI sparse arbitrary-depth EB composite-integral gates

Evaluate the complete physical four-level branching tree directly from sparse
owner fields after owner-local chemistry. Require roundoff-level parity with
the complete serial composite integral and require the global contributing-node
count to equal five.

Repeat for every valid level/patch selector and require parity with the serial
recursive subtree reference plus the exact topology-derived descendant count.
On multiple ranks, give rank zero a different but valid sibling selector; on
one rank, use an invalid selector. Require collective rejection, a zero
integral, and zero public local-node accounting. Retain the full established
MPI gates and all 208 serial tests in GNU Fortran Release and bounds/FPE-
checked Debug.

## 0.164.0 MPI owner-local arbitrary-depth EB hydro gates

Advance the post-chemistry physical four-level branching tree through one
small stable hydro interval after owner rotation. Require global per-level
node advances `[1, 2, 8, 8]`, and require grouped sender traffic to equal the
topology/owner-map sum of `parent_invocations * (refinement_ratio + 4)` for
every distinct-owner edge.

Materialize only after the transaction and require state/temperature parity
with the complete serial recursive hydro reference within qualified roundoff.
Also compare sparse and serial composite conserved integrals. On multiple
ranks, supply a different but valid interval on rank zero; on one rank, use a
negative interval. Require collective preflight rejection, exact sparse
rollback, zero advances, and zero transfers. Retain the full established MPI
gates and all 208 serial tests in GNU Fortran Release and bounds/FPE-checked
Debug.

## 0.165.0 MPI owner-local arbitrary-depth EB transport gates

Advance the post-hydro physical four-level branching tree through one small
stable SSPRK2 transport interval after owner rotation. Require global per-level
Euler advances `[2, 4, 16, 16]`. Require grouped sender traffic to equal twice
the qualified hydro edge-route formula plus one final distinct-owner
restriction transfer per relation.

Materialize only after the transaction and require state, temperature, and
minimum-limiter parity with the complete serial recursive transport reference
within qualified roundoff. Compare sparse and serial composite conserved
integrals independently. On multiple ranks, supply a different but valid
interval on rank zero; on one rank, use a negative interval. Require collective
preflight rejection, exact sparse rollback, unit limiter, zero advances, and
zero transfers. Retain all established MPI gates and all 208 serial tests in
GNU Fortran Release and bounds/FPE-checked Debug.

## 0.166.0 MPI owner-local arbitrary-depth EB full-physics gates

Starting from the qualified post-transport physical tree, apply one small
`R-T-H-T-R` interval with chemistry and every explicit-transport component
enabled. Require global chemistry advances `[2, 2, 4, 2]`, transport Euler
advances `[4, 8, 32, 32]`, and hydro advances `[1, 2, 8, 8]`. Require two
chemistry restriction schedules, two complete SSPRK2 transport schedules, and
one hydro schedule from the topology/owner-map traffic formulas.

Materialize only after the outer transaction. Require state, temperature,
minimum-limiter, and composite-integral parity with the complete serial split
reference within qualified accumulated roundoff. A rank-dependent or negative
interval must reject before mutation with unit limiter and zero counters in
every category. Retain all established MPI gates and all 208 serial tests in
GNU Fortran Release and bounds/FPE-checked Debug.

## 0.167.0 MPI owner-local arbitrary-depth EB clock gates

Start the qualified post-full-physics sparse and serial trees at time zero and
request a target equal to half the preceding full-physics test interval. Require
both clocks to take one exactly clipped step, publish the target time, report
the same minimum dt, and evaluate exactly five global owner nodes for timestep
selection.

Require the full-physics advance/transfer schedules from `0.166.0`, field and
temperature parity within `6e-9`, limiter parity within `1024 epsilon`, and
composite-integral parity within `6e-9`. Give rank zero a different valid target
time on multiple ranks and a negative target on one rank; require preflight
rejection with unchanged state/clock and neutral diagnostics. Repeat with a
zero step ceiling and require the same rollback. Retain all MPI gates and all
208 serial tests in GNU Fortran Release and bounds/FPE-checked Debug.

## 0.168.0 serial arbitrary-depth EB tagged-rebuild gates

Start from an eight-by-eight root-only EB tree containing one thermodynamically
consistent hot cell. Require normalized temperature tags to produce one
deterministic child at each of two relations, with both children attached to
their actual parent and the plan stopping at the requested three-level ceiling.

Apply the public tagged-regrid transaction and require a valid `1/1/1` tree,
roundoff-level composite-integral conservation, and a field-exact no-op when
the same plan is requested again. Reject the geometry builder and require the
accepted tree to remain exact. Finally make every node uniform, require zero
tags and a transactional collapse to one root level, and recheck the composite
integral. Retain all 208 serial tests in GNU Fortran Release and bounds/FPE-
checked Debug.

## 0.169.0 MPI owner-local arbitrary-depth EB tagged-rebuild gates

Start from the same eight-by-eight root-only hot-cell tree used by the serial
gate. Require the serial and owner-local planners to produce the same two
relations, parent indices, rectangles, and EB geometry while exactly one owner
evaluates every prospective parent.

Apply the public sparse regrid and require a valid `1/1/1` tree, exact
topology-derived remote prolongation/restriction counts, serial field and
temperature parity, and roundoff-level composite conservation at one, two,
four, and eight ranks. Expand the tag criteria to change retained overlap and
repeat every parity check. The unchanged plan must then be a field-exact no-op.

Give one rank different valid criteria on multiple ranks and invalid criteria
on one rank; require preflight rejection with exact sparse rollback and neutral
outputs. Finally make the hierarchy uniform, require a tag-free collapse to
the root, and recheck serial parity and conservation. Retain every established
MPI gate and all 208 serial tests in GNU Fortran Release and bounds/FPE-checked
Debug.

## 0.170.0 serial arbitrary-depth EB checkpoint gates

Write the qualified four-level branching tree after changed-topology overlap
retention. Read it into an empty candidate and require the same ordered
topology, lifecycle metadata, conserved fields, and recovered temperatures
within formatted roundoff.

Read the same file with a level ceiling below its stored depth and with the
first two species exchanged. Both operations must reject with an empty tree
and neutral metadata. Supply inconsistent time/step metadata to the writer and
require rejection before the valid checkpoint is replaced. Retain all 208
serial tests in GNU Fortran Release and bounds/FPE-checked Debug.

## 0.171.0 sparse MPI arbitrary-depth EB checkpoint gates

Write the qualified physical four-level branching tree through a selected I/O
root. Require the global sender count to equal the number of numerical nodes
not owned by that root. Read on the root, request a different work exponent,
recompute ownership, scatter directly, and require the sender count to equal
the number of nodes newly owned away from the root.

Materialize only after restart and require topology, state, recovered
temperature, time, minimum timestep, step count, and regrid count parity with
the serial source. Give one rank a different maximum depth on multiple ranks;
on one rank, request less depth than stored. Also exchange species only on rank
zero. Require both cases to reject with empty distribution/tree, zero metadata,
and zero transfers. Run at one, two, four, and eight ranks in GNU Fortran
Release and bounds/FPE-checked Debug, then retain all MPI gates and all 208
serial tests.

## 0.172.0 arbitrary-depth EB composite output gates

Write the qualified four-level branching serial tree to one CSV. Derive the
expected composite size from all node cell counts minus every child coarse-
rectangle area. Require the header to expose EB metrics and ordered species,
and require exactly one data row per expected leaf cell.

Repeat from sparse owner-local fields with the last rank as writer root.
Require the global sender count to equal the number of nodes not owned by that
root, and let only that root inspect and delete the file. Run this gate at one,
two, four, and eight ranks in GNU Fortran Release and bounds/FPE-checked Debug,
then retain all established MPI gates and all 208 serial tests.

## 0.173.0 runnable serial arbitrary-depth EB application gates

Run the public patch-tree executable from the established reactive 2D, EB, and
AMR namelists. Start from a 12-by-12 plane-EB hotspot and require recursive
temperature tags to populate levels zero through three before advancing to the
requested final time.

Read the single composite CSV and require unique level/patch/local-index cell
identities, the exact refinement-scaled spacing on every level, finite fields,
valid volume fractions, positive density/pressure/temperature, all three EB
cell classes, seven ordered species columns, and mass-fraction closure within
`8e-12`. Retain the complete 210-test suite in GNU Fortran Release and bounds/
FPE-checked Debug.

## 0.174.0 public patch-tree restart gates

Run the same four-level dynamically tagged case as an uninterrupted reference,
a process stopped after its first scheduled checkpoint, and a separate restart
process. Require the stopped composite time to be strictly inside the requested
interval and both completed paths to reach the exact final time.

Validate the checkpoint magic, schema, seven-species header, four-level count,
and completion marker. Key both completed composite outputs by
`(level, patch, i, j)`, require identical key and column sets, and compare every
numeric field within `3e-10 * max(1, abs(reference))`. Retain all 214 tests in
GNU Fortran Release and bounds/FPE-checked Debug.

## 0.175.0 public sparse-MPI application gates

Run the public arbitrary-depth EB input case at one, two, four, and eight MPI
ranks. Require four populated levels, identical identity-keyed composite
topology and columns, finite positive thermodynamics, seven-species closure,
the exact final time, and every numeric field within
`3e-10 * max(1, abs(one-rank reference))`. Retain the complete MPI gate chain
and all 214 serial tests in GNU Fortran Release and bounds/FPE-checked Debug.

## 0.176.0 public sparse-MPI cross-rank restart gates

Run an uninterrupted one-rank four-level application as the reference. In a
separate process, run the same physical case on two ranks with uniform node
weighting, stop after its first scheduled checkpoint, and require its
composite time to lie strictly inside the requested interval. Restart that
single checkpoint independently on four and eight ranks with depth-squared
weighting and require both continuations to reach the exact final time.

Validate the checkpoint magic, schema, seven-species header, four-level count,
and completion marker. Key the reference and both restarted composites by
`(level, patch, i, j)`, require identical key and column sets, and compare
every numeric field within `3e-10 * max(1, abs(one-rank reference))`. Retain
the complete MPI gate chain and all 214 serial tests in GNU Fortran Release
and bounds/FPE-checked Debug.

## 0.177.0 owner-local public sparse-MPI startup gates

For every fresh public application process, require the distribution to exist
before numerical initialization and require exactly one rank to execute the
root reactive initializer. The sparse root initializer must reject numerical
input allocated on any non-owner, validate the complete owner field, transfer
both allocatable arrays into the sparse node without a copy, and leave both
source arrays unallocated on every rank.

Exercise that path through the established one-, two-, four-, and eight-rank
four-level application parity gate. Retain the independent two-rank
checkpoint-stop and four-/eight-rank restart gate, the complete MPI gate
chain, and all 214 serial tests in GNU Fortran Release and bounds/FPE-checked
Debug.

## 0.178.0 public patch-tree checkpoint fingerprint gates

Write the public four-level checkpoint with schema 2 and require the structured
mesh/EB/physics/regrid fingerprint. Restart it through the existing serial and
MPI inputs and retain exact final topology plus complete-field parity. Change
only CFL in a separate restart input and require serial and four-rank MPI
processes to reject it before publishing topology, fields, or output.

Continue to permit the established rank-count and work-exponent changes across
the MPI checkpoint boundary. Retain the schema-1 low-level checkpoint tests,
the complete MPI gate chain, and all 215 serial tests in GNU Fortran Release
and bounds/FPE-checked Debug.

## 0.179.0 interface-local multilevel EB closure gates

Construct the recipient mask independently from child rectangles and parent EB
geometry. Require every selected cell to be active, unrefined, and within the
clipped three-by-three support of a direct coarse/fine interface; require no
selected cell inside any sibling rectangle.

Run the existing cut-interface fixed-depth, multipatch, branching patch-tree,
transport, chemistry, and public application gates. Require unchanged
composite mass, total energy, species closure, positive recovered temperature,
rollback behavior, and serial/MPI field parity in all 215 serial tests and the
complete one-, two-, four-, and eight-rank MPI Debug/Release chain.

## 0.180.0 embedded-wall transport gates

Recover one general-EOS H2/O2 state beside a unit-normal wall. Require a hotter
isothermal wall to add positive energy with exactly zero mass, momentum, and
species transfer. Require a moving no-slip wall to apply signed normal,
tangential, and out-of-plane momentum transfer plus positive wall work. Switch
the same wall to slip and require an exact zero viscous flux; supply a negative
centroid distance and require transactional rejection.

Compare the complete EB transport right-hand side with adiabatic-slip and
isothermal-no-slip wall records. Require differences only in EB cut cells and
nonzero energy and momentum changes there. Retain all single-level,
fixed-depth, multipatch, arbitrary-depth, serial, and sparse-MPI transport,
rollback, and field-parity gates in all 215 serial tests and the complete
one-, two-, four-, and eight-rank MPI Debug/Release chain.

## 0.181.0 public embedded-wall input gates

Run the public single-level EB transport application from a namelist selecting
an isothermal, tangentially moving no-slip wall. Require every active output
state to remain finite with species closure, require at least one cut cell to
heat relative to the no-transport reference, and require nonzero tangential
velocity in the cut-cell band. Retain the existing bulk-conduction reduction
of the temperature span.

Through the direct driver API, select isothermal mode with transport disabled
and require neutral clock/counter rejection. Keep active nondefault walls
rejected by checkpoint-capable AMR application preflight. Retain all 215 serial
tests and the complete one-, two-, four-, and eight-rank MPI Debug/Release
chain.

## 0.182.0 AMR wall-control restart gates

The fixed-depth restart formats compare the wall kind, thermal mode,
temperature, velocity, transport enable flag, each molecular-transport process
flag, and transport CFL before accepting state. The arbitrary-depth serial and
sparse-MPI applications carry the same values in their shared fingerprint.

Parity evidence requires an active isothermal wall to advance through a public
two-level AMR lifecycle, write and restart its checkpoint, and reject a changed
wall temperature without advancing the clock. The installed application also
reads a moving no-slip wall from its namelist and must produce cut-cell heating
and tangential momentum in the established two-level transport regression.
Exact PeleC checkpoint-file compatibility is not claimed; the contract is
PeleF-internal and intentionally rejects older schemas.

## 0.183.0 EB-safe limited-linear prolongation gates

Populate the coarse level with an EOS-consistent conserved state multiplied by
a nonconstant linear function of both Cartesian indices. Require an interior
regular child to reproduce the analytic MC-limited value and require
average-down of the complete fine patch to recover every parent state within
roundoff.

For every cut parent in the same geometry, require all children to retain the
qualified PCM value. Supply a nonfinite coarse component and require failure
with exactly zero state and temperature outputs. Run the complete serial suite
in GNU Fortran Release and bounds/FPE-checked Debug before qualifying the
kernel; retain the MPI chain before merging the versioned milestone.

## 0.184.0 fixed-depth public prolongation-selection gates

Set `prolongation_method = "linear"` in the installed hot-wall AMR transport
case and require its existing active-cell, thermal-span, moving-wall momentum,
and conservation checks to pass. Through the direct public driver, run linear
initialization and regridding in two-level, separated sibling-patch, and
three-level lifecycles.

Supply an unknown method to the shared dispatcher and require neutral outputs.
Attempt linear fixed-depth checkpointing and a linear arbitrary-depth
fingerprint and require rejection before time, counters, topology, or fields
are published. Retain all 215 serial regressions and the complete one-, two-,
four-, and eight-rank MPI Debug/Release chain.

## 0.185.0 restart-safe arbitrary-depth prolongation gates

Run every public fixed-depth checkpoint family with linear prolongation:
root/single-patch, separated sibling patches, static three levels, and dynamic
three levels. Require uninterrupted and checkpoint/restart outputs to retain
their existing field parity. Change only the method when reading a direct
two-level checkpoint and require neutral rejection.

Run the four-level public patch-tree case with linear prolongation in serial
and with one, two, four, and eight MPI ranks. Write a schema-4 checkpoint with
two ranks and restart it independently with four and eight ranks under a
different ownership weighting; require identity-keyed topology and field
parity against the uninterrupted one-rank reference. Change only the method in
the incompatible input and require rejection before topology or fields are
published. Retain all serial and MPI Debug/Release gates.

## 0.186.0 conservative limited-linear cut-parent gates

Use a diagonal-plane EB geometry and an EOS-consistent conserved field with
variation in both Cartesian directions. Require at least one cut parent to
produce nonconstant active fine children. Bound every active cut-parent child,
component by component, by the active 3-by-3 coarse-neighbor envelope.

Average the complete fine patch back with the established EB volume weights
and require every covered parent state to recover within roundoff, including
cut parents. Retain analytic regular-parent linear reproduction, temperature
recovery, nonfinite-input neutral rejection, and the parent-local PCM retry.
Run the complete serial and MPI Debug/Release chains before merging the
versioned milestone.

## 0.187.0 multidimensional cut-parent prolongation gates

Populate every active coarse state with one conserved vector multiplied by an
interface-tangential affine function evaluated at its fluid-volume centroid.
Prolong through the diagonal-plane EB patch and require every active child of
every cut parent to match the same analytic function at its fine fluid
centroid within roundoff.

Retain the nonconstant cut-parent, active component-envelope, complete-patch
average-down, regular-parent analytic value, EOS recovery, nonfinite rollback,
and parent-local PCM retry gates. Run the complete serial and one-, two-,
four-, and eight-rank MPI Debug/Release chains before merging.

## 0.188.0 rank-recovering cut-parent stencil gates

Construct a valid coarse/fine EB pair whose cut parent has only a straight
face-connected fluid path in its 3-by-3 neighborhood. Turn that path inside
the surrounding 5-by-5 box and populate every coarse cell with an EOS-valid
two-dimensional affine conserved field. The compact fit is rank one, while
the grown connected fit is full rank.

Require every active fine child to reproduce the analytic affine value,
including a nonzero difference in the direction unresolved by the compact
stencil. Average the child block down and require exact recovery of the parent
state. Retain the diagonal-plane, envelope, EOS/PCM retry, nonfinite rollback,
serial, and one-, two-, four-, and eight-rank MPI Debug/Release gates.

## 0.189.0 transactional fixed three-level parent-regrid gates

Initialize an offset reactive hotspot on a fixed three-level EB hierarchy, then
plan a different root-to-middle rectangle from the root temperature field.
Require deepest-first restriction, changed parent bounds, valid nested patch
descriptors, and the established two-cell finest safety margin after rebuild.

Compare every component of the three-level composite conserved integral before
and after the transaction. Supply an invalid prolongation control and require
neutral rejection with all three state/temperature fields and both patch
descriptors unchanged. Retain the complete serial and one-, two-, four-, and
eight-rank MPI Debug/Release chains before merging.

## 0.190.0 public dynamic-parent lifecycle gates

Enable the opt-in parent policy in a public offset-hotspot run and require the
accepted root-to-middle descriptor to differ from its configured seed. Require
a valid nested finest descriptor, a positive regrid count, complete EB cell
coverage, active thermodynamic validity, species closure, and composite
conservation.

Write a schema-4 dynamic checkpoint after a moved parent. Require the stored
policy bit, actual parent and finest descriptors, cadence counter, and terminal
marker. Restart in a separate process, derive expected output sizes from the
stored descriptors, and compare every root, middle, and finest field against
an uninterrupted run. Retain legacy fixed-parent inputs and the complete
serial and MPI Debug/Release chains.

## 0.191.0 arbitrary-depth outflow-boundary child gates

Run a public four-level x-upper hotspot in the serial patch-tree application.
The composite checker requires every populated level to reach the exact
physical boundary, refinement-ratio-consistent cell spacing, all EB classes,
finite positive thermodynamics, and species closure.

Run the same input through the sparse-MPI application at one, two, four, and
eight ranks. Check the one-rank topology contract and require every ordered
composite identity and numeric field to match across rank counts before the
ordinary Debug/Release and complete MPI regression chains are accepted.

## 0.192.0 boundary-touching patch-tree restart gates

Write a schema-4 four-level x-upper checkpoint after at least one committed
step. Require the stopped composite output and every independent continuation
to retain exact boundary contact at all levels.

Compare the serial uninterrupted and restarted outputs by composite identity
and every numeric field. Repeat with a two-rank checkpoint and four- and
eight-rank restarts under a different ownership weight. Retain terminal-marker,
fingerprint-mismatch, topology, time, and field-parity gates in the complete
Debug/Release chains.

## 0.193.0 boundary-touching recursive transport gates

Enable Fourier conduction in the public fresh and split-run x-upper trees.
Require four populated levels, exact boundary contact, finite positive fields,
species closure, and serial/sparse 1/2/4/8-rank composite parity after the
recursive transport/hydro transaction.

Require the schema-4 checkpoint to retain the active transport fingerprint.
Compare independent serial continuation and two-rank to four-/eight-rank
sparse continuation against their uninterrupted references, then retain the
complete Debug/Release regression chains.

## 0.194.0 boundary-touching mixture-transport gates

Run the existing x-upper fresh and split-run trees with viscosity, Fourier
conduction, mixture-averaged species diffusion, and barodiffusion all enabled.
Retain four populated levels, exact physical-side contact, finite positive
fields, species closure, and serial/sparse 1/2/4/8-rank composite parity.

Require the same controls in the schema-4 checkpoint fingerprint. Compare the
independent serial continuation and two-rank to four-/eight-rank sparse
continuations against uninterrupted full-transport references, including the
changed ownership weight.

## 0.195.0 reacting boundary-tree full-physics gates

Enable elementary chemistry beside all four transport controls in the public
x-upper fresh and split-run inputs. Require the complete transactional
`R-T-H-T-R` clock to retain four populated levels, exact boundary contact,
finite positive thermodynamic fields, species closure, and serial/sparse
1/2/4/8-rank composite parity.

Require schema 4 to retain the chemistry model and tolerances with the
transport fingerprint. Compare independent serial continuation and two-rank
to four-/eight-rank sparse continuation against uninterrupted reacting
references under the changed ownership weight.

## 0.196.0 checkpointed transport-limiter gates

Write a nonneutral minimum transport theta through the direct serial and
selected-root sparse checkpoint APIs. Require exact round-trip recovery,
schema-5 public headers, communicator-wide agreement, and neutral value `1`
on rejected reads.

Run the public reacting x-upper checkpoint chain unchanged. Its serial and
two-to-four/eight-rank continuations must retain uninterrupted fields while
the application resumes its accumulated limiter minimum from checkpoint.

## 0.197.0 checkpointed conservation-baseline gates

Write a synthetic finite `nvar`-component baseline through the direct serial
and selected-root sparse APIs. Require exact schema-6 round-trip recovery,
collective rejection when ranks supply different values, and no allocated
baseline after a rejected read.

Capture the uninterrupted and changed-rank restart application diagnostics.
In addition to identity-keyed field parity, require the restarted final
conservation error and minimum transport limiter to match the uninterrupted
logical run.

## 0.198.0 checkpointed AMR operator-counter gates

Write synthetic chemistry, transport, and hydro vectors whose capacity is
larger than the stored topology depth. Require exact schema-7 serial and
selected-root sparse recovery, rejection when the restart depth cannot hold
the vector, collective rejection when ranks disagree, and no allocated
counter output after a rejected read.

For public applications, capture the three reported level-counter lines from
the uninterrupted run and every independent restart. Require exact vector
equality after the serial split and after two-rank checkpoint continuation on
four and eight ranks, in addition to the existing field, conservation, and
limiter gates.

## 0.199.0 checkpointed AMR regrid-history gates

Write synthetic regrid-evaluation and cumulative tagged-cell values through
the direct serial and selected-root sparse APIs. Require exact schema-8
recovery, rejection of partial or inconsistent history, communicator-wide
rejection of rank disagreement, and neutral scalar outputs after failed reads.

Capture both public diagnostics in the existing serial and sparse restart
logs. Require the independent serial continuation and two-rank checkpoint to
four-/eight-rank continuations to match their uninterrupted reference exactly.

## 0.200.0 public branching patch-tree gates

Use separated boundary-touching and interior hotspots to create multiple
ordered patches on at least one populated level. Inspect composite identities
to reject a degenerate one-patch chain while retaining four populated levels,
x-upper contact at every level, finite positive fields, and species closure.

Run the same reacting full-transport case in serial and sparse MPI at 1/2/4/8
ranks. Reuse the serial and two-to-four/eight-rank restart chains and require
identity-keyed topology and every numeric field to match the uninterrupted
branching reference. The schema remains 8.

## 0.202.0 three-dimensional directional-core gates

Verify exact z-to-x-to-z round trips for conserved and primitive states.
Compare z-normal equal-state Rusanov and PeleC-style fluxes with the analytic
Euler physical flux, including both transverse momenta. Require passive-species
z fluxes to sum to the mass flux and reacting H2/O2 z-normal HLLC fluxes to
preserve analytic momentum, energy, and species placement. Independently
verify x/y/z mesh spacings and endpoint cell centers.

These are local algebraic and thermodynamic gates. Whole-step 3D dimensional
reduction and convergence remain mandatory before the regular 3D solver can be
called qualified.

## 0.203.0 conservative periodic 3D Euler gates

Initialize the same smooth constant-pressure state as x-, y-, and z-dependent
fields. Rotate the y- and z-dependent results back to canonical x layout and
require the complete PCM/SSPRK2 update to match the established 1D periodic
solver for both Rusanov and PeleC-style fluxes. Reject an invalid solver
without changing any caller-owned state.

Advect a diagonal periodic entropy wave on 16^3, 32^3, and 64^3 grids. Require
monotonically decreasing density L1 error, observed orders above 0.85, positive
density and pressure, and relative conservation error below `2e-11` for mass,
all three momenta, and total energy. Run the public 16^3 input separately and
check its row count, coordinate cardinality, exact-wave L1 bound, primitive to
conserved relations, and analytical integral totals from the emitted CSV.

These are analytical and cross-dimensional internal gates, not a claim of
whole-application parity with a three-dimensional PeleC executable.

## 0.204.0 general-EOS multispecies 3D Euler gates

Build one elementary H2/O2/N2 general-EOS line state and advance an independent
periodic 1D SSPRK2 reference directly from the public x-normal flux. Embed and
rotate that line through x, y, and z uniform 3D fields. Require the complete
state and recovered temperature to agree for Rusanov, HLLC, and PeleC-style
fluxes. Reject an invalid solver without changing either caller-owned array.

Advect a constant-composition, constant-pressure diagonal entropy wave on
6^3, 12^3, and 24^3 grids. Require density L1 order above 0.70, positive
density/pressure/temperature, species closure below `5e-11`, and relative
conservation below `5e-11` for every Euler and species component. Separately
run the public 12^3 full-H2O2-thermo application and check its dynamic species
schema, x-fastest ordering, exact-wave L1 bound, ideal-mixture EOS,
primitive/conserved identities, composition, and analytical mass, momentum,
and species totals.

These gates qualify composition-dependent 3D hydrodynamics only. They do not
exercise reaction rates, chemistry splitting, or molecular transport and do
not establish whole-application PeleC parity.

## 0.205.0 cell-local 3D chemistry gates

Initialize uniform elementary and full-H2O2 3D fields and advance one chemistry
interval. Require every cell's conserved state and recovered temperature to
match the established one-dimensional constant-volume reactor path. Require
nonzero species evolution and H/O/N elemental conservation for both the
explicit seven-species and implicit ten-species mechanisms.

For a uniform field, require the complete 3D Strang transaction to reduce to
two independent one-dimensional chemistry half-steps because periodic uniform
hydro has zero divergence. Let the first reaction half-step succeed and then
select an invalid hydro solver; the outer transaction must retain every
caller-owned state and temperature value exactly.

Run the public periodic 6^3 elementary-H2/O2 Gaussian hotspot. Independently
check its dynamic species schema, row and coordinate order, physical fields,
EOS and conserved/primitive identities, nonnegative mass fractions, closure,
analytical initial mass, H/O/N elemental totals, zero net momentum, nonzero
chemistry evolution, and resolved temperature structure. These are internal
coupling and conservation gates, not external ignition validation or
whole-application PeleC parity. Molecular transport remains unqualified.

## 0.206.0 regular-grid 3D molecular-transport gates

Embed the same nonuniform reacting line in x-, y-, and z-dependent uniform
3D fields. Advance one full transport SSPRK2 interval with viscosity,
conduction, species diffusion, and barodiffusion enabled. Rotate z momentum
back to the canonical line layout and require elementary seven-species and
full-H2O2 ten-species state, temperature, and limiter results to match the
established 2D operator to component-scaled roundoff. Independently verify
the summed `1/dx^2 + 1/dy^2 + 1/dz^2` stable-step formula and exact rollback
for an incompatible transport table.

Diffuse a periodic diagonal three-dimensional divergence-free shear mode on
6^3, 12^3, and 24^3 grids. Require both observed L1 orders above 1.70,
inactive species limiting, and roundoff total-energy and momentum
conservation. Separately require Fourier conduction to contract a temperature
span and mixture-averaged/barodiffusive species transport to reduce H2
mass-fraction variance while conserving every periodic species integral and
retaining nonnegative species density.

Require transport-disabled full composition to match the established
`R-H-R` result exactly. With transport enabled, reject an invalid middle hydro
step after a valid transport prefix without changing caller state or
temperature. Run paired public 6^3 full-H2O2 inert hotspots and require the
transport result to retain EOS/closure identities and every conserved
integral while producing resolved thermal and species differences from the
control.

These gates qualify only the selected serial, single-level, periodic transport
subset. They are not evidence for Soret/Dufour or multicomponent transport,
AMR reflux, restart, MPI decomposition, EB walls, or whole-PeleC parity.

## 0.207.0 static two-level reactive 3D AMR gates

Require PCM prolongation followed by restriction to reproduce every covered
parent exactly, composite integration to avoid double-counting covered coarse
cells, and average-down to change no uncovered cell. Require the new
flux-returning coarse SSPRK2 interface to reproduce the established update
exactly and to zero its output without changing caller fields on rejection.

Advance a genuinely three-dimensional elementary-H2/O2 entropy wave through a
strictly interior ratio-two patch. Require a nonzero reflux correction,
roundoff conservation of every composite Euler and species integral, exact
covered-parent synchronization, positive density/pressure/temperature,
species closure, and complete rollback for an invalid face solver. A separate
uniform-field case must remain invariant across subcycling and synchronization.

Run the public full-H2O2 case through `pelef_amr_reactive_3d`. Independently
parse both level CSV files, validate dynamic species relationships and
physicality, and reconstruct every covered coarse conserved component from
the `2^3` fine children. These gates do not qualify chemistry or transport on
AMR, high-order reconstruction, dynamic regridding, restart, MPI rank parity,
boundary-touching refinement, EB, or whole-application PeleC parity.

## 0.208.0 static 3D AMR checkpoint/restart gates

Round-trip a synchronized elementary-H2/O2 hierarchy through the public
checkpoint API. Require exact conserved arrays and scalar metadata, recovered
temperature agreement, and rejection of unsynchronized writes. Change the
solver fingerprint and truncate the file independently; each read must fail
without changing any caller array, temperature, time, step, integral, or
reflux value. Configuration tests reject checkpoint or restart paths that
would overwrite a CSV output.

For the public full-H2O2 ratio-two case, run an uninterrupted reference, stop
after coarse step four, and resume from the resulting checkpoint to the same
final time. Compare the complete coarse and fine CSV byte streams, not only
norms or selected fields. Exact identity is required. These gates establish
serial same-build continuation for the static single-patch hydro hierarchy;
they do not establish schema migration, cross-toolchain bitwise identity,
distributed restart, dynamic patch restoration, or physics beyond the
qualified AMR boundary.

## 0.209.0 distributed-slab static 3D AMR gates

Partition both coarse and fine levels into contiguous x slabs and run the
kernel with one, two, four, and eight OpenMPI ranks. Each rank must own at
least one plane at the eight-rank edge. Require the distributed CFL value,
one complete AMR advance, both temperature fields, and the maximum reflux
correction to match the serial kernel bit for bit.

Collective preflight must reject a rank-disagreed root, field extent, solver,
patch descriptor, floating control, or NASA7 species value before entering a
variable-count collective. Every such rejection must leave state,
temperature, continuation metadata, and reflux output unchanged.

Run the public full-H2O2 case continuously in serial and at 1/2/4/8 ranks.
Apply the independent physical/coarse-fine checker to the MPI result and
require byte-for-byte equality of both level CSV streams at every rank count.
Then write a synchronized checkpoint at coarse step four with two ranks and
resume it independently with four and eight ranks. Both continuations must
match the uninterrupted serial files byte for byte.

These gates qualify rank-owned arithmetic over a replicated hierarchy and
selected-root formatted I/O. They do not establish sparse persistent fields,
distributed checkpoint layout, performance scaling, AMR source/diffusion
coupling, dynamic topology, high-order reconstruction, 3D EB, or complete
PeleC parity.

## 0.210.0 regular and static-AMR characteristic-PLM gates

Rotate one nonuniform general-EOS line through x, y, and z fields and require
the complete PLM SSPRK2 state and recovered temperature to agree after the
corresponding momentum rotation. Unknown reconstruction, limiter, or Riemann
controls must leave state and temperature unchanged and return zero face
fluxes.

Advect the elementary-H2/O2 diagonal entropy wave on 6^3, 12^3, and 24^3
periodic grids. Require positive state, species closure, all-component
conservation, and observed density orders above 1.65. Run a separate
full-H2O2 namelist case through the public executable and independently check
its x-fastest CSV state relationships and analytical density error.

Advance the static hierarchy with PCM and characteristic PLM. Require
roundoff composite conservation, exact average-down synchronization, physical
coarse/fine states, nonzero reflux, and a PLM fine-grid density L1 no more than
70 percent of the PCM error. Exercise both MC and minmod on a uniform
hierarchy. A step-four PLM checkpoint must resume to coarse and fine CSV files
that are byte-identical to an uninterrupted run; changing reconstruction in
the restart target must reject transactionally.

These gates qualify a serial method-of-lines characteristic PLM path and
limited coarse ghost prolongation. They do not establish CTU/PPM, distributed
PLM halo exchange, physical boundaries, AMR source/diffusion coupling,
dynamic topology, or 3D EB.

## 0.211.0 static axis-plane 3D EB gates

Construct x-, y-, and z-normal planes analytically and require exact fluid
volume, cut-cell count, cut fraction, cell and face centroids, EB area,
physical centroid, and solid-to-fluid normal. In every cell, independently
reconstruct the signed Cartesian aperture divergence and compare it with the
stored EB normal integral. Mutating either side of that identity must make the
geometry invalid. Invalid axes, bounds, extents, nonfinite positions, and an
interior grid-face-aligned plane must be rejected.

Recover an elementary-H2/O2 mixture pressure from conserved state and require
the 3D slip wall to emit only the negative pressure traction along the
solid-to-fluid normal. Its result must be independent of velocity. Integrate
the wall source over x/y/z plane cuts and require the expected unit-area
pressure force with no mass, energy, or species leakage. Cartesian equal-state
pressure fluxes plus the wall source must leave every active cell at zero
right-hand side to roundoff, while covered cells stay exactly zero.

Advance uniform flow tangent to the x/y/z planes using HLLC, Rusanov, and
PeleC-style PCM faces. Require exact state and temperature invariance. Advance
a localized density perturbation with Rusanov and require roundoff
fluid-volume conservation of every component, species-density closure,
positive temperature, and unchanged covered cells. Unknown solvers, invalid
timesteps, nonphysical cut states, nonfinite fluxes, and bad extents must fail
transactionally.

These gates establish an internal analytic geometry and conservation
reference; they are not a direct field comparison to the pinned PeleC build.
They do not establish oblique/curved geometry, face-aligned EB ownership,
small-cell redistribution, high-order cut-face reconstruction, chemistry or
transport at EB cells, AMR/reflux/regridding, restart, MPI, or a public 3D EB
application.

## 0.212.0 conservative planar 3D EB FluxRedist gates

Apply arbitrary residuals to fully regular and fully covered geometry and
require exact identity and zero respectively. For x-, y-, and z-normal planes
with `kappa=0.05`, compare the cut and one receiving-cell values against the
closed-form FluxRedist result. Require less than eleven percent of the raw cut
residual, exact uniform active-residual preservation, compact support, and
fluid-volume conservation for every component.

On an x-normal `8 x 3 x 3` cut sheet, seed the central cut cell and require
nonzero transfer through the positive-x face and both y and z directions.
The overlapping neighborhoods must retain the original weighted integral.
Invalid geometry and nonfinite residuals must return a zero unpublished
result.

Advance an elementary-H2/O2 state from a supplied residual that would make
the raw cut state negative. Require positive cut and receiving states after
redistribution, unchanged temperature for a uniformly scaled state, exact
component conservation, and transaction rollback for a larger nonphysical
request.

Finally construct the complete PCM Rusanov residual for a `kappa=0.05`
density jump. At a three-dimensional CFL-0.5 full-grid step, independently
prove that the raw cut density is negative and that the raw update rejects.
Require the redistributed path to remain physical and conserve all Euler and
species integrals to roundoff. Repeat x/y/z tangent-flow invariance and all
invalid solver/timestep/state transactions through the integrated path.

These are internal formula, conservation, and stability gates against the
frozen 2D FluxRedist semantics. They are not a direct PeleC field comparison
and do not qualify weighted StateRedist, arbitrary geometry, an EB-aware
public timestep selector, high-order reconstruction, coupled physics, AMR,
restart, MPI, or an application-level 3D EB workflow.

## 0.213.0 zeroth-order planar 3D EB StateRedist gates

For x-, y-, and z-normal `kappa=0.05` planes, set one provisional cut state to
`-U` in an otherwise uniform active field. With target `0.5`, independently
derive neighbor weight `0.45`, cut factor `0.6363636363636364`, and receiver
factor `0.9181818181818182`; require all components to match. Require compact
support, uniform-state identity, covered standalone zero, and
fluid-volume-weighted conservation.

The geometry selector must follow the sole nonzero aperture-difference-normal
component to a regular receiver. A domain-edge small cell without that
receiver, an invalid target, a nonfinite provisional state, or any unsupported
normal must fail without publishing a partial result.

From a caller-supplied reactive residual, require StateRedist to repair a
provisional negative cut density, recover its NASA7 temperature, and conserve
all Euler/species components. Setting `target=kappa` must remove the merge and
therefore expose the nonphysical provisional state; a larger bad update must
also roll back state and temperature exactly.

Finally feed raw, FluxRedist, and StateRedist from one common PCM face/wall
divergence. Repeat all x/y/z uniform tangent-flow gates, nonuniform
conservation/species closure/covered identity, and invalid
solver/timestep/state/target transactions. At the frozen CFL-0.5 full-grid
control, independently require raw negative density and rollback while both
named stabilization routes remain physical and conservative.

These are analytical and invariant-based internal gates for order-zero
StateRedist on exact axis planes. They do not establish higher-order or
general-geometry AMReX/PeleC StateRedist parity, a public 3D EB solver,
coupled physics, AMR/restart, or MPI ownership.

## 0.214.0 public stabilized planar 3D EB gates

Recover the frozen elementary-mixture primitive state on every active cell and
compare the EB CFL result with the direct x/y/z spectral-rate formula. Covered
storage must not affect the result. Invalid CFL, an inadmissible active state,
bad extents, and an all-covered geometry must return zero timestep.

Validate all three plane orientations and an actual namelist read. Reject
nonfinite configuration, face alignment, a cut cell without a regular
receiver, invalid StateRedist target, unsupported solver/method, and public
selection of the raw route.

Run the installed app twice on the same `10 x 8 x 6`, `kappa=0.05`,
constant-pressure density sheet, selecting StateRedist and FluxRedist. Both
must take two positive steps to `5.0e-5`, remain physical, preserve mass,
energy, species, and tangential momenta to roundoff, and report the normal
wall-momentum exchange separately.

Independently parse both 480-row CSV files. Reconstruct x-fastest indices,
Cartesian centers, 144/48/288 covered/cut/regular classifications, fractions,
and fluid centroids. Recompute mixture-EOS and conserved/primitive relations,
analytical mass/species totals, resolved density/pressure response, and
cross-method invariant agreement. Freeze one distinct maximum-density and
normal-impulse signature for each algorithm.

These gates establish a reproducible public workflow, not external PeleC
field parity or physical validation. They do not cover higher-order EB hydro,
coupled chemistry/transport, general geometry, AMR/restart, or MPI execution.

## 0.215.0 optional elementary chemistry in planar 3D EB gates

Advance x-, y-, and z-oriented 3D fields through the cell-local elementary
reactor with an active mask. Require roundoff conservation of density, total
energy, and H/O/N totals, a nonzero species response, exact identity of one
masked regular cell and every covered cell, and whole-field rollback for an
invalid mask shape.

Compose masked chemistry around the planar EB hydro transaction. The direct
unit gate must contain a `kappa=0.05` cut cell below the `0.5` StateRedist
target, preserve mass, total energy, both tangential momenta, and H/O/N totals,
change at least one active species, and retain covered state and temperature
bitwise. Chemistry-disabled execution must match the direct StateRedist result
bitwise. Empty chemistry, invalid solver, and invalid reactor tolerances must
reject without publishing any candidate.

Run public inert and reacting StateRedist cases beside the frozen `0.214.0`
StateRedist baseline. Require the inert CSV to be byte-identical to the
baseline. Independently parse the reacting 480-row CSV, reconstruct geometry,
EOS and species closure, verify all 144 covered rows column-for-column against
the inert output, preserve mass/energy/tangential momentum and H/O/N totals,
and freeze both the active-species-change and maximum-density signatures.

These gates qualify only the serial, single-level, first-order, pre-reaction-
CFL-controlled case with the fixed seven-species elementary mechanism. They do
not establish general reactive timestep stability, a production stiff
mechanism, public FluxRedist chemistry, molecular transport, higher-order or
general-geometry EB, AMR/restart, MPI ownership, or external PeleC parity.

## 0.216.0 molecular transport in planar 3D EB gates

For x-, y-, and z-normal `kappa=0.05` planes, exclude covered storage from the
transport-step scan and require the same positive bound after overwriting that
storage with invalid values. Evaluate nonuniform velocity and composition
fields through the aperture-weighted transport residual and SSPRK2 StateRedist
advance. Require exact-zero outer and embedded-wall face fluxes, a resolved
state response, all-component and H/O/N conservation to roundoff, exact
covered identity, and activation of the species-inventory limiter under a
forced oversized interval.

The full split unit must prove that transport-disabled execution is bitwise
identical to the `0.215.0` `R-H-R` wrapper. With transport and chemistry active,
require a resolved coupled response, fluid and elemental conservation, exact
covered storage, and a bounded limiter diagnostic. Invalid solver,
redistribution, process dependency, transport-table extent, interval, or
post-stage EOS state must reject without publishing any part of the candidate.

Run public control, transport-only, and coupled StateRedist cases on the same
`10 x 8 x 6`, `kappa=0.05`, constant-pressure density sheet. Stop after one
clipped `2.0e-9 s` step so the zero-gradient hydro outer boundary does not
enter the conservation comparison. Parse all 480 rows and require 144 covered
rows to remain byte-identical. Independently reconstruct geometry, NASA7 EOS,
energy, closure, physical-volume component and H/O/N integrals, active
transport/chemistry response, operator flags and sequence, step schedule,
transport diagnostics, and the complete Debug/Release CSV SHA-256 set.

The frozen signatures are transport response `1.0489647406832604e-1`,
chemistry response `4.043407673691575e-8`, effective relative conservation
error `6.867e-15`, elemental error `1.112e-14`, maximum diffusivity
`9.7180419287635654e-4`, and limiter factor `1.0`. These are internal formula,
transaction, and invariant gates, not direct PeleC field parity or external
validation. They do not qualify nonzero physical-boundary transport, general
embedded-wall models, FluxRedist transport, higher-order/general geometry,
AMR/restart, MPI, or long-time reacting-flow stability.

## 0.217.0 transactional planar 3D EB restart gates

Round-trip a physical coupled planar EB state through the schema-1 checkpoint
and require exact state, temperature, inventories, clock, and cumulative
diagnostic recovery. Mutate immutable configuration, reaction, and transport
records independently. Truncate the file and make a stored temperature
inconsistent with its conserved state. Every negative read must reject and
leave all caller-owned arrays and scalar metadata bitwise unchanged.

Run the installed application continuously and as two processes separated
after the first committed chemistry-plus-transport step. The public case must
take at least one additional step after restart. Require byte-identical final
CSV output, identical final step/time and cumulative minimum-step/diffusivity/
limiter diagnostics, an earlier distinct checkpoint-stop CSV, and explicit
magic, schema, and terminal-marker records.

These tests establish deterministic serial continuation under the exact stored
physics contract. They are not an external PeleC field comparison, a scalable
checkpoint performance result, a general schema-compatibility guarantee, or
evidence for MPI/AMR/general-geometry restart.

## 0.218.0 pinned Cantera ingestion gates

Hash the repository YAML and require explicit selection of `ohmech`, not the
co-located Redlich-Kwong phase. Import it twice and require deterministic JSON.
Verify ten ordered species, 29 ordered reactions, NASA7 at 101325 Pa, complete
gas transport, element balance, eleven third-body records, one Troe record,
six duplicate flags, and source indices 1 through 29. Negative fixtures reject
NASA9, missing transport, PLOG, Chebyshev, SRI, Lindemann-only falloff,
unbalanced reactions, invalid identifiers, nonfinite values, and oversized
species sets.

Regenerate the committed normalized JSON under Cantera 3.2 and compare it
byte-for-byte. Then regenerate the Fortran module with the dependency-free
Python generator and compare it byte-for-byte in every build configuration.
The compiled full-H2/O2 support test checks the embedded source provenance,
reaction-family counts, duplicate/order metadata, thermo/transport ordering,
and a positive mixture-transport evaluation.

Finally rerun the existing Cantera production-rate/0D trajectory gate and the
regular 1D/2D full-H2/O2 parity pair. Those comparisons are the numerical
guard against a semantically wrong named-collider or unit conversion. They do
not establish arbitrary-mechanism parsing, detailed-fuel validation, runtime
selection, or CVODE parity.

## 0.219.0 selected mechanism build gates

Generate a Fortran module from a normalized two-species fixture as an ordinary
build dependency, compile it into an isolated library, and run a configured
probe. The probe must load ordered NASA7, reaction, and primitive transport
records, select a temperature inside their common range, and obtain finite
production rates and a finite nonzero Jacobian. The fixture deliberately uses
a third-body reaction so this is more than a syntax-only module import.

For a user-selected bundle, configure a clean tests-disabled Release build
with the pinned full-H2/O2 JSON and its YAML source. Require configure-time and
build-time source hash checks, a successful ten-species/29-reaction probe, an
installed probe byte-identical to the build-tree executable, and a
non-executable stack. Pair that with a negative configuration using the wrong
YAML and require rejection before generation.

Treat the selected bundle and source as configure dependencies. Touch copied
inputs and require a plain build to rerun CMake, recheck the source, regenerate,
and relink. Compile a maximum-length 128-character reaction equation with all
generated lines bounded to 100 characters, and reject invalid or colliding
Fortran names before publication.

Compare the energy-constrained reduced reactor Jacobian with a centered
directional finite difference of the full thermochemical RHS. This prepares a
stronger callback gate for a future stiff backend without claiming CVODE
integration in this milestone.

These gates qualify one configure-time normalized-bundle ABI. They do not
establish runtime mechanism dispatch, CFD coupling of the selected bundle,
new thermo/rate families, mechanisms above 32 species, or physical validation.

## 0.223.0 selected regular 1D reacting-flow gates

Build the selected two-species fixture through the same CMake mechanism-bundle
helper used by an external bundle. Run a uniform reacting case with chemistry
and molecular transport enabled. Independently parse its dynamic species CSV
schema and require the configured bundle order, finite positive thermodynamic
state, mass-fraction closure, spatial uniformity, nonzero chemical activity,
the requested final time, and the application's zero reported invariant drift.
Repeat the run and require a byte-identical file.

Build the pinned full-H2/O2 bundle as a separate selected runtime. Run the
ordinary fixed and configured selected regular-1D applications from matching
inputs and require byte-for-byte CSV identity, in addition to independent
schema, physicality, closure, final-time, and uniformity checks. This is a
stronger gate than tolerance comparison because it exercises the same
hydrodynamic, native reaction, and mixture-transport implementation through
two separately linked mechanism-loading paths.

Give the selected application an unknown composition species and give the
fixed application `chemistry_model = "selected"`; both must fail during
startup and create no output artifact. Unit gates also reject duplicate names,
nonfinite unused composition entries, and a composition-wave endpoint that
would become negative, while proving safe normalization of finite near-maximum
fractions.

Finally rerun the eight Debug/Release, Cantera, MPI, and optional SUNDIALS
configurations. In a clean tests-disabled selected Release build, run all three
selected executables, install every application, require build/install byte
identity and non-executable GNU stacks, and verify that the installed selected
1D executable has no RPATH/RUNPATH or unresolved library. These gates establish
one configure-time serial regular-grid 1D dispatch boundary. They do not
establish runtime loading, selected 2D/3D/AMR/EB/MPI coupling, CVODE in CFD,
PeleC field parity, detailed-fuel validity, or production performance.

## 0.224.0 selected regular 3D reacting-flow gates

Build the selected two-species fixture through the public mechanism-bundle
helper and run a 4-by-4-by-4 periodic uniform reactor with chemistry and
molecular transport enabled. Independently parse its dynamic species CSV and
require bundle ordering, finite positive density/pressure/temperature,
mass-fraction closure, spatial uniformity, nonzero H2 conversion, the requested
final time, and zero reported mass, momentum, and energy drift.
Repeat the run and require a byte-identical file.

Build the pinned full-H2/O2 bundle in a separate selected runtime. Run matching
fixed and selected periodic regular-3D transport-hotspot inputs and require
byte-for-byte CSV identity. Retain the existing independent checker requiring
positive state, composition closure, conservation, and a nonzero H2 transport
response relative to the uniform control. This separates loader/configuration
parity from the shared time-integration driver without claiming two independent
3D algorithms.

Run a short nonuniform ten-species/29-reaction full-H2/O2 hotspot with
characteristic PLM and chemistry enabled through both fixed and selected front
ends. Require measurable H2 conversion, positive state, closure, the requested
time, and byte-exact fixed/selected CSV identity. This is a dispatch and active
reaction-family smoke gate, not an ignition-delay or physical-validation case.

Give the selected front end an unknown composition species and give the fixed
front end `thermo_model = "selected"`; both must fail during startup without a
CSV artifact. Unit gates also exercise duplicate-name, active/unused-tail,
finite-value, and overflow-safe normalization behavior through the shared
composition resolver. Generator tests render and compile the selected 3D
template using maximum-length legal generated identifiers.

Finally rerun the eight Debug/Release, Cantera, MPI, and optional SUNDIALS
configurations. In a fresh tests-disabled selected Release build, run all four
selected executables, install every application, require build/install byte
identity and non-executable GNU stacks, and verify that the installed selected
3D executable has no RPATH/RUNPATH or unresolved library. These gates establish
configure-time serial periodic single-level regular-3D dispatch. They do not
establish physical-boundary, selected 2D/AMR/EB/MPI, runtime loading, CVODE in
CFD, PeleC field parity, detailed-fuel validity, or production performance.

## 0.225.0 selected regular 2D reacting-flow gates

Build the two-species selected fixture through the public bundle helper and
run a 4-by-4 uniform 2D reactor with chemistry and transport enabled.
Independently parse the dynamic species schema and require bundle ordering,
positive density/pressure/temperature, closure, spatial uniformity, measurable
H2 conversion, and the requested final time. Repeat the run and require
byte-identical output.

Run matching fixed and selected pinned full-H2/O2 2D inputs in three regimes.
The uniform `2.0e-6 s` case must exercise the generic adaptive implicit
chemistry path and produce byte-identical CSV. A `4 x 4`, `1.0e-9 s`
nonuniform hotspot must enable characteristic PLM, CTU, and all 29 reactions;
require positive state, closure, measurable H2 conversion, and byte-exact
fixed/selected output. An `8 x 8`, `5.0e-7 s` physical-boundary case must
enable species diffusion and a zero-net-mass prescribed wall flux and again
produce byte-identical output.

The selected parser must remain opt-in. Reject selected input sent to the
fixed executable and reject an unknown bundle species before output creation.
Unit gates cover duplicate names, nonfinite composition tails, overflow-safe
normalization, selected composition-wave endpoints, a valid actual-species
wall flux, and a nonzero wall-flux slot beyond the actual bundle size. Runtime
negative gates require conservative double-hotspot and isothermal reflected
ghost-temperature ranges to remain within the common NASA7 interval.

Generator tests must render and compile the selected 2D template using
maximum-length legal identifiers. Rerun the eight standard configurations and
perform a fresh tests-disabled selected Release build/install audit with all
five selected executables. Require build/install byte identity, non-executable
stacks, no selected-2D RPATH/RUNPATH or unresolved libraries, and installed
smoke/parity runs. These gates establish configure-time serial single-level
regular-2D dispatch, including the existing physical-boundary set. They do not
establish selected AMR/EB/MPI, runtime loading, CVODE inside CFD, external
PeleC field parity, detailed-fuel validity, or production performance.

## 0.226.0 selected serial AMR 1D reacting-flow gates

Build the two-species selected fixture through the public bundle helper and
run a nonuniform reactive hotspot with chemistry, molecular transport, and
solution-driven ratio-two refinement enabled. Independently parse the dynamic
species schema and require levels zero and one, positive state, mass-fraction
closure, exact requested time, level-width ratio, gap-free/nonoverlapping
composite domain coverage, measurable H2 conversion, and at least one fine
region. Repeat the run and require byte-identical output.

Run matching fixed and selected pinned full-H2/O2 inputs in three regimes. A
two-level nonuniform hotspot must enable chemistry and molecular transport and
produce byte-identical CSV. An arbitrary-depth case must create three levels
with characteristic PPM and preserve exact fixed/selected output. A dynamic
multipatch transport case must create at least three disjoint fine regions and
again produce byte-identical output.

The selected parser must remain opt-in. Reject an unknown bundle species,
non-AMR input, an entropy wave whose constant-pressure density amplitude sends
`T0/(1-A)` outside the common NASA7 range, and selected input sent to the fixed
AMR executable before output creation. Generator tests must reserve the new
program/module namespace and render and compile the selected AMR template with
maximum-length legal identifiers.

Rerun the eight standard configurations and perform a fresh tests-disabled
selected Release build/install audit with all six selected executables.
Require build/install byte identity, non-executable stacks, no selected-AMR
RPATH/RUNPATH or unresolved libraries, and installed fixture/full-H2/O2
two-level, multilevel, and multipatch smoke/parity runs. These gates establish
configure-time serial reactive-1D AMR dispatch. They do not establish selected
AMR checkpoint/restart, EB/MPI, runtime loading, CVODE inside CFD, external
PeleC field parity, detailed-fuel validity, or production performance.

## 0.227.0 selected serial reactive EB 2D gates

Build the two-species selected fixture through the public bundle helper and
run a plane-cut 4-by-4 case with native chemistry and molecular transport.
Independently parse the dynamic species schema and require exactly 8 regular,
4 cut, and 4 covered cells, valid volume fractions, positive active state,
mass-fraction closure, exact requested time, measurable H2 conversion, and
repeat-byte output identity.

Run matching fixed and selected pinned full-H2/O2 inputs in two regimes. A
plane-cut case must execute active chemistry and produce byte-identical CSV. A
circle-cut characteristic-PLM case must execute viscosity, conduction,
species diffusion, barodiffusion, correction velocity, species enthalpy flux,
and an isothermal moving no-slip embedded wall; require nonzero H2 spatial
response and byte-identical fixed/selected CSV.

The selected parser must remain opt-in. Reject an unknown bundle species, an
embedded-wall temperature outside the common NASA7 interval, a lexical
input/output alias, and selected input sent to the fixed EB executable before
output creation. Unit gates must require fixed-parser rejection and exact
caller-output rollback after an invalid selected chemistry-integrator policy.
Generator tests must reserve the new namespace and render and compile the
selected EB template with maximum-length legal identifiers.

Rerun the eight standard configurations and perform a fresh tests-disabled
selected Release build/install audit with all seven selected executables.
Require build/install byte identity, non-executable stacks, no selected-EB
RPATH/RUNPATH or unresolved libraries, and installed full-H2/O2 chemistry and
transport smoke/parity runs. These gates establish configure-time serial
single-level reactive 2D EB dispatch. They do not establish selected EB AMR/3D,
checkpoint/restart, MPI, runtime loading, CVODE inside CFD, external PeleC
field parity, detailed-fuel validity, or production performance.

## 0.228.0 selected regular MPI 1D gates

Extract the fixed 19-cell periodic decomposition, deterministic initialization,
adaptive `R-T-H-T-R` loop, collective invariants, ordered gather, and CSV
publication into one mechanism-independent driver. Freeze the old fixed
one-rank CSV before extraction and require the refactored fixed executable to
remain byte-identical on 1/2/4 ranks.

Build independently ordered two-species and pinned full-H2/O2 selected bundles
through the public mechanism helper. The two-species case must use explicit
chemistry, complete to `1.0e-7 s`, change both temperature and a species
integral, and produce byte-identical 1/2/4-rank CSV. The full bundle must use
implicit chemistry and remain byte-identical across ranks and to the fixed
MPI output. Compare dynamic `rhoY_*` columns rather than assuming H2/O2 names.

First run a nonzero valid adaptive-Strang control and require a chemical state
change. Then run an invalid policy at one rank or valid explicit/implicit
policies that disagree across ranks on 1/2/4 ranks; require collective failure
plus exact caller state and temperature rollback. Reserve the new Fortran
namespaces, render and MPI-compile the selected template with maximum legal
identifiers, and audit the selected link graph for absence of `pelef_core` and
committed mechanisms.

Rerun the eight standard configurations and perform a fresh tests-disabled
selected MPI Release build/install audit with all 36 installed executables.
Require build/install byte identity, non-executable stacks, no selected-MPI
RPATH/RUNPATH or unresolved libraries, a single coherent MPI implementation,
and installed fixed/selected 1/2/4-rank parity. These gates do not establish a
general input-driven MPI application, selected MPI AMR/EB, restart, runtime
loading, CVODE inside CFD, thread safety, scaling, detailed-fuel validity, or
external validation.

## 0.229.0 selected serial EB AMR 2D gates

Extract a shared fixed/selected application driver without changing the
existing static two-level arithmetic. Run the established fixed EB AMR
regression after extraction. Require the selected runtime link graph to omit
`pelef_core`, fixed database loaders, and committed generated mechanisms.

Run an independently ordered two-species static plane fixture twice. On the
coarse level require 32 regular, eight cut, and 24 covered cells; on the fine
level require 84 regular, 12 cut, and 48 covered cells. Require finite positive
state, mass-fraction closure, the exact final time, H2 change above `1.0e-5`,
exact ratio-two coordinates, volume-fraction consistency, composite
conserved-state average-down below `5.0e-12`, and byte-identical repeated CSV
independently on both levels.

Run matching fixed and selected pinned full-H2/O2 hotspot cases with
characteristic PLM, native implicit chemistry, viscosity, thermal conduction,
species diffusion, barodiffusion, and an isothermal moving no-slip embedded
wall. Require nonzero H2 change, coarse/fine average-down consistency, and
byte-identical fixed/selected coarse and fine CSV. Reject unknown species, an
out-of-range reflected wall temperature, every input/coarse/fine path alias,
dynamic, three-level, multipatch, and dynamic-parent modes, each checkpoint
control, restart, and selected input sent to the fixed executable. Preserve
the input and create neither level output on every startup rejection.

Reserve the new Fortran namespaces and render and compile the selected
template with maximum-length legal identifiers. Rerun the eight standard
configurations and perform a fresh tests-disabled selected MPI Release
build/install audit with all 37 installed executables. Require build/install
byte identity, non-executable stacks, no RPATH/RUNPATH, coherent MPI runtime
dependencies, and installed fixed/selected full-H2/O2 EB AMR parity. These
gates do not establish selected dynamic/three-level/multipatch EB AMR,
restart, selected EB 3D, MPI EB AMR, runtime loading, CFD CVODE, thread safety,
scaling, detailed-fuel validity, or external validation.

## 0.230.0 selected serial EB 3D gates

Freeze the fixed coupled seven-species Debug and Release CSV before extracting
one mechanism-independent application driver. Require the refactored fixed
wrapper to retain those files byte-for-byte. Build a complete test-only
selected bundle from the same reaction, thermo, and transport records and
require its selected output to be byte-identical to the fixed output in both
configurations.

Run an independently ordered H2/H bundle twice on a 4-by-4-by-4 exact
non-face-aligned plane. Require 32 regular, 16 cut, and 16 covered cells,
finite positive active state, mass-fraction closure, active H2 change, exact
covered storage, and byte-identical repeated CSV. Run the pinned full-H2/O2
bundle on a 10-by-8-by-6 plane with implicit chemistry and every qualified
molecular-transport term active; require 288 regular, 48 cut, and 144 covered
cells, nonzero H2 and density changes, positivity, and closure.

Reject unknown or duplicate species, nonfinite/zero/negative composition,
out-of-range initial temperatures, lexical input/output aliases, checkpoint
controls, restart, unsupported diagnostic element families, and selected
input sent to the fixed executable before output creation. Reserve the new
Fortran namespaces, compile the maximum-length rendered template, and audit
the selected target graph for absence of `pelef_core` and committed
mechanisms.

Rerun the eight standard configurations and perform a fresh tests-disabled
Release build/install audit with all 38 executables. Require build/install
byte identity, non-executable stacks, no RPATH/RUNPATH, resolved dependencies,
one coherent MPI ABI, and installed fixed and selected EB 3D smoke/checks.
These gates do not establish selected EB3D checkpoint/restart, arbitrary
element diagnostics, general geometry, EB AMR/MPI 3D, runtime loading, CFD
CVODE, thread safety, scaling, detailed-fuel validity, or external validation.

## 0.231.0 selected serial EB 3D restart gates

Freeze the Debug and Release schema-1 checkpoint-stop artifacts before changing
the adapter. After adding selected persistence, rerun the same fixed executable
and require the checkpoint and stopped CSV to retain their build-type-specific
exact SHA-256 values in the automated checker. The first common-body marker
must remain `SPECIES`; no selected record may be written by a fixed caller.

At the checkpoint API, write and read complete selected schema 2. Require the
`SELECTED_CONTEXT` marker, canonical bundle SHA-256, exact generated integrator,
species count, normalized bundle-order composition, and the unchanged complete
model/geometry/state body. Reject schema 1 through a selected reader, schema 2
through a fixed reader, incomplete context, changed SHA, changed integrator,
changed composition, malformed context records, negative/nonfinite composition,
truncation, wrong target shape, and existing model/configuration/state
incompatibilities. Every failed read must leave all caller geometry, arrays,
clock, counters, inventories, and cumulative diagnostics unchanged.

Run the complete seven-species selected elementary executable in three separate
processes on the established 24-by-8-by-6 coupled planar case: uninterrupted,
stop after the first committed checkpoint, and restart to completion. Require
the stopped output to precede final time, two committed final steps, exact
uninterrupted/restarted 1152-row CSV bytes, exact cumulative summaries, and a
schema-2 checkpoint containing the generated bundle identity and explicit
policy. Require the selected uninterrupted reference to equal the matched
fixed reference byte-for-byte. A fourth process must reject a changed normalized composition without
creating output. Reject lexical input/checkpoint and input/restart aliases
before model loading or file replacement.

Retain the fixed schema-1 restart chain and the complete eight-configuration
matrix, then repeat the tests-disabled selected Release build/install and ELF
audit. These gates establish provenance and transactional compatibility, not a
cryptographic digest of the checkpoint payload, crash-atomic replacement,
schema migration, scalable I/O, general EB geometry, AMR/MPI EB 3D restart,
physical validity of a selected mechanism, or external PeleC field parity.

## 0.232.0 selected static two-level EB AMR 2D restart gates

Freeze the fixed schema-3 checkpoint and checkpoint-stop coarse CSV before
adding selected context. Require their exact SHA-256 values in the fixed
restart checker after the change. Fixed calls must continue to write the same
header, first species record, patch, coarse/fine body, and terminal marker.

At the two-level API, require selected schema 4 to contain
`SELECTED_CONTEXT`, the canonical full-bundle SHA-256, `implicit` policy,
ten-species count, and normalized bundle-order composition before the
unchanged body. Exercise selected round trip, partial read/write context,
both cross-schema directions, changed SHA, policy, and composition, malformed
marker, truncation, invalid composition, failed-write non-destruction, and
complete rollback. Require failure text to distinguish missing context,
cross-schema input, and each selected identity mismatch.

Run five separate full-H2/O2 processes on the same static 8-by-8 coarse and
12-by-12 fine plane hierarchy: fixed uninterrupted, selected uninterrupted,
selected first-step checkpoint-stop, selected restart, and changed-composition
restart. Require 32/8/24 coarse and 84/12/48 fine regular/cut/covered counts,
two final coarse steps, a one-step checkpoint, exact fixed/selected reference
bytes, exact uninterrupted/restarted bytes on both levels, and no output from
the mismatch process. Freeze the selected reference coarse/fine, stopped
coarse/fine, and checkpoint SHA-256 values. Reject input/checkpoint and
restart/output lexical aliases before any output is created.

Rerun the full configuration matrix and repeat the tests-disabled selected
Release build/install and ELF audit. These gates establish deterministic
serial static two-level continuation and selected-context provenance only;
they do not qualify dynamic, three-level, multipatch, or MPI EB AMR
persistence, payload authentication, crash-atomic/scalable I/O, detailed-fuel
physics, or external PeleC parity.

## 0.233.0 selected dynamic two-level EB AMR 2D restart gates

Rerun the frozen fixed schema-3 and selected static schema-4 restart chains
after introducing dynamic selected persistence. Their exact checkpoint and
CSV hashes must not change. At the API level, require a selected dynamic call
to write and read schema 5, retain the complete selected context and initial
composite-integral baseline, and reject schema 4; require a selected static
call to reject schema 5. Corrupt the baseline marker and require rejection.
Every rejection must preserve caller-owned patch, arrays, clock, counters, and
diagnostics.

Use a 12-by-12 full-H2/O2 coarse grid whose initial `(2:5,2:5)` fine patch
does not contain the temperature hotspot. Regrid once to `(5:12,4:12)` before
the first-step checkpoint, producing a 16-by-18 fine grid. Run fixed and
selected uninterrupted references, selected checkpoint-stop, independent
selected restart, and changed-composition restart. Require three final coarse
steps, one checkpointed step, one completed regrid, exact fixed/selected
reference bytes, exact uninterrupted/restarted bytes on both levels, and no
output from the mismatch process.

The checker must validate schema 5, selected context, species order, dynamic
flags, stored post-regrid patch, regrid metadata, row counts, cell-type counts,
physicality, mass-fraction closure, logs, and frozen reference/stopped/
checkpoint SHA-256 values. Require uninterrupted/restarted maximum composite
conservation-error identity. Rerun the complete configuration matrix and a fresh
tests-disabled Release install with installed fixed and selected dynamic
restart smoke. These gates do not qualify selected three-level, multipatch,
dynamic-parent, or MPI EB AMR persistence, payload authentication,
crash-atomic/scalable I/O, detailed-fuel physics, or external PeleC parity.

## 0.234.0 selected static three-level EB AMR 2D restart gates

Rerun both frozen fixed three-level restart chains after introducing selected
static persistence: static magic/schema 3 and dynamic magic/schema 4 hashes
must not change. At the API level, require a selected static three-level call
to write and read schema 4 under the static magic, retain the complete selected
context and initial composite-integral baseline, and reject fixed schema 3.
The fixed reader must reject selected schema 4. Corrupt the baseline marker,
change bundle SHA, integrator, or composition, and require every failure to
publish no candidate and return topology, arrays, clock, counters,
diagnostics, and baseline in their empty/default failure states.

Run fixed and selected full-H2/O2 uninterrupted references, a selected
first-step checkpoint-stop, independent selected restart, and a changed-
composition restart on an 8-by-8 root with 12-by-12 middle and 16-by-16 finest
levels. Require two final coarse steps, one checkpointed step, exact fixed/
selected reference bytes, exact uninterrupted/restarted bytes on all three
levels, and no output from the mismatch process.

The checker must validate selected context, baseline, species order, topology,
time/step metadata, 64/144/256 row counts, cell-type counts, physicality,
mass-fraction closure, log identity, and frozen reference/stopped/checkpoint
SHA-256 values. Require uninterrupted/restarted cumulative conservation text
identity. Also reject lexical aliases involving the finest output before any
write. Rerun the complete configuration matrix and a clean tests-disabled
Release install with installed fixed and selected restart smoke. These gates
do not qualify selected dynamic-three-level, multipatch, dynamic-parent, or
MPI EB AMR persistence, payload authentication, crash-atomic/scalable I/O,
detailed-fuel physics, or external PeleC parity.

## 0.235.0 selected dynamic three-level EB AMR 2D restart gates

Rerun both frozen fixed three-level chains and the selected static chain.
Fixed dynamic magic/schema 4 and all earlier checkpoint and CSV hashes must
remain unchanged. At the API level, require a selected dynamic three-level
call to write and read schema 5 under the dynamic magic, retain the complete
selected context and original composite-integral baseline, and reject schema
4. Require the fixed dynamic reader to reject schema 5, and require both
static/dynamic selected magic directions to fail without publishing a
candidate.

Corrupt the baseline marker and numeric baseline, dynamic controls, committed
middle patch, first state field, and terminal marker. Reject an invalid
baseline write without replacing an existing checkpoint. Every failure must
leave topology, arrays, clock, counters, diagnostics, and baseline in their
documented empty/default states.

Run fixed and selected full-H2/O2 uninterrupted references, a selected
first-step checkpoint-stop, independent selected restart, and a changed-
composition restart on a 12-by-12 root. Require dynamic-parent regridding to
move the middle patch from `(2:11,2:11)` to `(5:10,3:10)` and the finest patch
from `(6:9,6:9)` to `(3:10,3:14)` before the checkpoint. The committed level
shapes must be 12-by-12, 12-by-16, and 16-by-24. Require two final coarse
steps, one checkpointed step, one completed regrid, exact fixed/selected and
uninterrupted/restarted bytes on all three levels, identical cumulative
conservation text, and no output from the mismatch process.

The checker must validate schema 5, selected context, composite baseline,
species order, dynamic controls, committed patches, clock/regrid metadata,
complete checkpoint fields, CSV topology and physicality, mass-fraction
closure, chemistry activity, logs, and frozen reference/stopped/checkpoint
SHA-256 values. Rerun the complete configuration matrix and a clean tests-
disabled Release install with installed fixed and selected dynamic-three-level
restart smoke. These gates do not qualify selected multipatch or MPI EB AMR
persistence, payload authentication, crash-atomic/scalable I/O, detailed-fuel
physics, or external PeleC parity.

## 0.236.0 selected dynamic multipatch EB AMR 2D restart gates

Rerun the frozen fixed patch-set restart chain and require its schema-3
checkpoint SHA-256 to remain
`3ce931f5dcd017de4864789f53b57a47e40fe82cc5fa6c87e9292ff7efe1022c`.
At the API level, require a selected patch-set call to write and read schema 4,
retain the complete selected context and original composite-integral baseline,
and reject schema 3. Require the fixed reader to reject selected schema 4.

Change composition, corrupt the baseline or terminal marker, append trailing
content, and reject an invalid baseline write without replacing an existing
checkpoint. Every read failure must leave the child-patch list, arrays, clock,
counters, diagnostics, and baseline in their documented empty/default states.
Pass an invalid selected integrator directly to the patch-set advance and
require rejection at the first root chemistry half-step, proving the policy is
not dropped before any child advance.

Run fixed and selected full-H2/O2 uninterrupted references, a selected first-
step checkpoint-stop, independent selected restart, and a changed-composition
restart on a 14-by-14 root with two disjoint 10-by-10 children. Require three
final coarse steps, one checkpointed step, one completed regrid, exact fixed/
selected reference bytes, exact uninterrupted/restarted bytes on root and both
children, identical cumulative conservation text, and no output from the
mismatch process.

The checker must validate schema 4, selected context, composite baseline,
species order, complete committed patch set, clock/regrid metadata, every
checkpoint field through end-of-stream, CSV topology and physicality, mass-
fraction closure, chemistry activity, logs, and frozen reference/stopped/
checkpoint SHA-256 values. Also require the selected front end to reject every
derived child-output alias before mechanism loading. Rerun the complete
configuration matrix and a clean tests-disabled Release install with installed
fixed and selected multipatch restart smoke. These gates do not qualify
selected MPI EB AMR persistence, payload authentication, crash-atomic/scalable
I/O, detailed-fuel physics, or external PeleC parity.

## 0.237.0 selected sparse MPI AMR 1D restart gates

Rerun the frozen fixed patch-tree regression and require its schema-1
checkpoint SHA-256 to remain
`1337b54d2761f2e3a2d25e264e9d6d073d1934756b206da16634d3140d933ed6`.
At the API level, require selected schema 2 to round-trip its complete context
and original composite baseline, reject schema 1, and leave caller state
unchanged after mismatch, truncation, or trailing content. Require the fixed
reader to reject schema 2 and an invalid selected write to preserve its target.

Run the sparse chemistry policy unit on one, two, and four ranks. An explicit
policy must be accepted; an invalid policy must fail collectively at the first
chemistry half-step with exact solution rollback. Give different ranks
different otherwise-valid compositions and require context consensus to abort
before output.

Run fixed one-rank and selected one-/two-/four-rank uninterrupted full-H2/O2
references, a selected one-rank first-step checkpoint-stop, and independent
two-/four-rank restarts. Require three active levels, four final coarse steps,
exact fixed/selected and rank-count bytes, exact uninterrupted/restarted bytes,
physicality, closure, active HO2/H2O2, and the frozen selected checkpoint and
CSV hashes. Mutate bundle, integrator, composition, species order, and baseline
independently. Also mutate density, raw species density, an EOS-tolerance-scale
negative raw species density, species closure, temperature, and geometry.
Require every restart to reject without input mutation or output. Reject
input/output aliasing before mechanism loading.

Rerun the complete configuration matrix and a clean tests-disabled Release
install with installed fixed and selected changed-rank restart smoke. These
gates do not qualify selected MPI EB AMR, payload authentication, crash-atomic
or scalable I/O, detailed-fuel physics, or external PeleC parity.

## 0.238.0 selected sparse MPI reactive EB patch-tree 2D gates

Preserve the fixed schema-8 sparse MPI EB lifecycle while moving it into one
shared fixed/selected driver. Under matched GNU Release flags, require the
new fixed checkpoint and stopped CSV to remain byte-identical to the final
`0.237.0` executable. Pin the Debug and Release checkpoint hashes separately
because floating-point optimization changes trailing digits.

Run a planar-EB sparse chemistry policy unit on one, two, and four ranks.
Reject an invalid one-rank policy and otherwise-valid rank-dependent policies
collectively before candidate mutation, then accept a uniform explicit policy.
Give different ranks different normalized selected compositions and require
the shared driver to abort before root initialization or output.

Run the fixed full-H2/O2 one-rank reference and selected one-, two-, and four-
rank processes on the same 12-by-12 dynamic case. Require levels 0 through 3,
active chemistry on every level, an initialization regrid and post-step regrid
evaluation, physical finite output, species closure, and byte-exact fixed/
selected/rank-count CSV parity. Require the selected front end to reject fixed
chemistry, input/output aliasing, checkpoint cadence/path, and restart path;
require the fixed front end to reject selected chemistry. No rejected process
may publish an output or checkpoint.

Also build the selected MPI EB target from the independent two-species fixture
bundle, run its generated explicit policy through the same four-level dynamic
lifecycle, and validate the resulting topology, EB classes, finite state, and
species closure independently of the full-H2/O2 exact-parity case.

Rerun the complete configuration matrix and a clean tests-disabled Release
install with installed selected one-/two-/four-rank smoke and ELF/MPI ABI
audit. These gates qualify configure-time selected execution only. They do not
qualify selected MPI EB persistence, fixed-depth MPI modes, runtime loading,
scalable I/O, physical detailed-fuel validity, or external field parity.

## 0.239.0 selected sparse MPI reactive EB patch-tree 2D restart gates

Keep the frozen fixed schema-8 Debug and Release checkpoint hashes. At the
serial checkpoint API, require schema-9 round trip with selected context and
original composite baseline, fixed-to-selected and selected-to-fixed rejection,
composition mismatch rejection, and invalid-baseline write non-replacement.

Run fixed one-rank and selected one-/two-/four-rank uninterrupted full-H2/O2
references. Stop a selected one-rank process after its first checkpoint and
resume independently on two and four ranks. Require four active levels, two
final root steps, exact fixed/selected/rank-count CSV bytes, exact
uninterrupted/restarted bytes, restored operator and regrid histories, finite
physical state, species closure, and frozen GNU Debug/Release schema-9,
stopped, and final hashes.

Mutate schema, species order, bundle, integrator, composition, numerical
fingerprint, geometry, composite baseline, density, raw species, near-floor
negative species, species closure, temperature, terminal marker, truncation,
and trailing content independently. Every restart must fail before output and
must leave its checkpoint bytes unchanged. Exercise all six lexical alias
classes among input, output, checkpoint, and restart paths before mechanism
loading. On two ranks, pass the selected fingerprint on only one rank at the
read and write boundaries; both calls must reject collectively before root I/O
or gather, transfer no entities, and publish no checkpoint.

Rerun the complete configuration matrix and a clean tests-disabled Release
install with the installed one-rank stop and two-/four-rank restart chain.
These gates establish rank-neutral selected persistence, not payload
authentication, crash-atomic/scalable I/O, schema migration, runtime loading,
fixed-depth MPI modes, performance, detailed-fuel physics, or external parity.

## 0.240.0 rank-local sparse static 3D AMR gates

At one, two, four, and eight ranks, require the parent-aligned distribution to
cover every coarse and fine x plane exactly. Check previous/next active fine
neighbors directly, require four zero-fine-cell ranks in the eight-rank case,
and prove that every multi-rank hierarchy stores fewer numerical values than
the complete global hierarchy. Compare both periodic coarse ghost layers with
the corresponding serial planes bit for bit.

Advance the coarse periodic level and the complete static hierarchy with
characteristic PLM. Require exact serial parity for conserved state,
temperature, all three time-averaged face-flux families, reflux magnitude,
and synchronized coarse/fine state. Repeat the public PCM and PLM entropy-wave
cases at one, two, four, and eight ranks and compare both CSV files byte for
byte with their serial references.

Reject invalid and otherwise-valid rank-dependent reconstruction, NASA7, root,
and patch/layout contracts before mutation or a variable-count collective.
Every failure must keep the sparse hierarchy bit-identical and zero the reflux
result. Write PCM and PLM checkpoints with two ranks, restart independently at
four and eight ranks, and require byte-exact serial final files. The schema
must contain neither communicator size nor ownership.

Run GNU Debug and Release focused sets, the complete configuration matrix, and
a tests-disabled installed MPI Release smoke. These gates qualify rank-local
static periodic inviscid hydro and rank-neutral compatibility I/O. They do not
qualify dynamic topology, AMR chemistry/transport, physical boundaries,
CTU/PPM, 3D EB AMR, scalable I/O, performance scaling, or external PeleC
field parity.

## 0.241.0 fixed static 3D AMR chemistry gates

At the serial API, advance an elementary uniform hierarchy through a source
phase and require species activity, exact Euler invariance, H/O/N conservation,
average-down synchronization, and rollback for invalid reaction, integrator,
and middle-hydro stages. Compare chemistry-disabled Strang dispatch directly
with the established hydro entry.

At one, two, four, and eight MPI ranks, require bitwise serial parity for the
source phase and the complete `R-H-R` split. Reject otherwise-valid
rank-dependent reaction records and chemistry-enable flags without mutation.
At eight ranks, corrupt coarse state only on a rank with no fine cells and
require collective source rejection with exact local rollback.

Run the full-H2/O2 hotspot and inert control through the public serial
application. Validate exact schema and row topology, state relationships,
physicality, closure, nonzero species and temperature response, Euler and
H/O/N conservation. Repeat the reacting application at one, two, four, and
eight ranks and compare both CSV levels byte for byte with serial. Run these
gates in GNU Debug and Release, retain all prior hydro/restart hashes, and
rerun the complete configuration matrix. These are structural, conservation,
and deterministic parity gates, not ignition or external field validation.

## 0.242.0 selected static 3D AMR gates

Run the generated two-species fixture twice in serial and require exact coarse
and fine bytes. Repeat it at one, two, and four MPI ranks, with the four-rank
layout proving that ranks owning no fine planes remain active. Apply the
generic output checker because H/O/N diagnostics are intentionally unavailable
for this mechanism.

Run generated full-H2/O2 through the selected serial frontend and at one, two,
four, and eight MPI ranks. Require both levels to match the fixed `0.241.0`
serial files byte for byte and require their frozen SHA-256 values. Retain the
fixed characteristic-PLM checkpoint, coarse-output, and fine-output hash gates
so a common fixed/selected regression cannot self-reference to green.

At two and four ranks, mutate one otherwise-valid selected context field on one
rank: SHA, integrator, composition, NASA7 coefficient, reaction equation,
reaction rate, or transport record. Every rank must return rejection without
abort or output. Public negative cases must reject selected transport,
checkpoint, restart, and lexical path aliases before creating output, and the
fixed frontend must reject selected input. Run focused Debug and Release gates,
the complete eight-configuration matrix, and a tests-disabled installed MPI
Release smoke. These gates do not qualify selected restart, AMR transport,
physical validation, performance, or scaling.

## 0.243.0 selected static 3D AMR restart gates

1. Preserve the fixed schema-2 checkpoint and characteristic-PLM checkpoint,
   coarse, and fine SHA-256 gates without routing fixed calls through schema 3.
2. Require all selected context fields before opening a checkpoint and reject
   fixed/selected cross-schema reads without changing caller state.
3. Round-trip schema 3 exactly and reject changed bundle SHA, integrator,
   normalized composition, reaction rate, or chemistry tolerance while state,
   temperatures, clock, step, baseline, and reflux history remain unchanged.
4. Reject malformed rank-local reaction efficiency shape collectively before
   shape-dependent packing at two and four ranks. Require chemistry tolerances
   to agree by representation before root I/O.
5. Stop a full-H2/O2 selected serial run after its first step, restart it, and
   compare both final CSV files byte-for-byte with an uninterrupted run.
6. Write the same checkpoint on two MPI ranks and resume it independently on
   one, two, four, and eight ranks. Require exact serial reference bytes and a
   zero-fine-plane diagnostic at eight ranks.
7. Require serial and MPI checkpoint byte identity and freeze schema-3,
   coarse, and fine hashes as `f44f69...c423`, `a03884...0ff3`, and
   `8d5a99...fcbe`.
8. Preserve input bytes and forbid outputs when checkpoint/input paths alias.

These gates prove deterministic rank-neutral continuation and provenance
binding. They do not authenticate payloads or qualify crash-atomic/scalable
I/O, AMR transport, dynamic topology, physical boundaries, performance,
ignition, detailed fuels, or external PeleC field parity.

## 0.244.0 static 3D AMR molecular-transport gates

1. Compare the hierarchy parabolic step with independent coarse and fine
   regular-grid limits and require `min(dt_c, r^2 dt_f)` exactly.
2. Advance a nonuniform serial hierarchy with all transport processes and
   require finite positive temperatures, species closure, average-down
   synchronization, nonzero reflux, and composite conservation. Repeat with
   viscosity, conduction, and species diffusion enabled individually and with
   refinement ratio three.
3. Force an extreme species gradient independently across x/y/z lower and
   upper coarse/fine faces. Require an active interface theta below one,
   nonnegative coarse and fine species, successful EOS recovery, and composite
   conservation after the conservative fine-side flux correction.
4. Reject negative, greater-than-one, and NaN external ghost-face theta values
   with zero outputs and unchanged inputs. Reject barodiffusion without species
   diffusion in configuration and numerical APIs.
5. Compose active chemistry and transport through `R-T-H-T-R`; require a
   synchronized physical result and exact rollback when the middle hydro stage
   rejects its candidate.
6. At one, two, four, and eight MPI ranks, compare transport timestep,
   diffusivity, theta, reflux, conserved fields, and recovered temperatures
   bit for bit with serial. Repeat the six active-interface cases; the
   eight-rank case must include ranks with no fine planes.
7. Reject a rank-local transport-record or transport-policy mismatch before
   mutation or payload-dependent communication. Require identical rollback on
   every rank.
8. Run fixed and selected full-H2/O2 public cases in serial and at one, two,
   four, and eight ranks. Check schema, physicality, species closure,
   conservation, resolved transport activity, synchronization, and frozen
   coarse/fine SHA-256 values.
9. Require transport-enabled checkpoint or restart input to fail before any
   output because schemas 2 and 3 do not persist the parabolic operator
   contract. Retain all prior fixed and selected restart hashes.

These gates qualify one periodic, static, strictly interior two-level
hierarchy and deterministic rank-local arithmetic. They do not qualify
dynamic topology, physical coarse boundaries, 3D EB AMR, scalable I/O,
performance, detailed-fuel physics, ignition, or external PeleC field parity.

## 0.245.0 static 3D AMR transport-checkpoint gates

1. Write and read fixed transport checkpoints as exclusive schema 4 for the
   public `full_h2o2` case with `chemistry_enabled=.false.` and
   `transport_enabled=.true.`. Write and read selected transport checkpoints
   as exclusive schema 5 with `thermo_model='selected'`,
   `chemistry_enabled=.true.`, and `transport_enabled=.true.`.
2. Require both schemas to bind the complete ordered gas-transport database,
   parameter convention, operator identity, all transport controls, and
   cumulative transport diagnostics. Require the phase to be
   `POST_ACCEPTED_COARSE_STEP` after a committed `R-T-H-T-R` step.
3. Validate private read candidates through the terminal marker and strict EOF
   before publication. Mutate each transport record/control/diagnostic and
   require rejection with caller state unchanged. Retain fixed schema-2 and
   selected schema-3 transport-disabled compatibility and reject cross-schema
   fallback.
4. Require all MPI ranks, including ranks with no fine planes, to agree on
   transport and selected context before payload communication or root I/O.
   Require root read/write failure to be reported collectively and prevent any
   rank from continuing.
5. Write a serial and a two-rank root-formatted checkpoint, resume it at the
   supported changed rank counts, and compare coarse/fine fields and restored
   transport diagnostics exactly with the uninterrupted reference.
6. The focused gates pass GNU Debug/Release serial `10/10` each and
   GNU/OpenMPI Debug/Release MPI `23/23` each. Freeze the schema-4/schema-5
   checkpoint and output hashes in
   [`validation/0.245.0.md`](validation/0.245.0.md).

These focused gates cover transport persistence for one static, strictly
interior, periodic, two-level Cartesian hierarchy. No 0.245 full configuration
matrix or clean installed Release result is claimed. Formatted replacement is
not crash-atomic; payload authentication, scalable I/O, dynamic topology,
physical boundaries, 3D EB AMR, performance, and external PeleC field parity
remain outside this increment.
