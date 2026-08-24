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
- equal-state physical/Rusanov/HLLC flux identity;
- stationary and moving heterogeneous-composition contact preservation;
- exact equality between summed species flux and total mass flux;
- homogeneous hydro update equal to zero;
- homogeneous Strang-split field equal to independent zero-dimensional cell chemistry;
- global mass, three momenta, and total-energy conservation;
- density, pressure, temperature, and species positivity;
- mass-fraction closure;
- smooth entropy-wave convergence above order 1.75 on both refinement intervals;
- smooth H2/N2 composition-wave convergence above order 1.70;
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
