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
