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
