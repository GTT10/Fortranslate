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

The current Cantera gate establishes parity only for four reversible elementary reactions without third-body or falloff effects. The reactive-flow tests establish numerical coupling and reduction properties for that same subset; they do not establish parity for Cantera's complete `h2o2.yaml`, a stiff mechanism, PelePhysics chemistry integration, molecular transport, or multidimensional PPM/transport parity. The periodic 2D CTU subset is qualified separately below.

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

The periodic general-EOS CTU path is accepted only while all of the following remain active:

- equal-state y-direction HLLC flux has the correct normal/tangential momentum placement;
- the sum of species fluxes equals the y mass flux;
- a y-uniform 2D step agrees with the 1D characteristic-PLM/HLLC update to roundoff;
- the oblique constant-pressure entropy wave converges on 12, 24, and 48 square grids;
- the transverse correction has a measurable signature and does not materially degrade the 24-square result;
- all corrected face states and final cells recover positive density, pressure, temperature, and closed composition through the NASA7 EOS;
- global mass, both in-plane momenta, out-of-plane momentum, and total energy remain conservative;
- a periodic velocity vortex remains positive and preserves its nontrivial velocity/pressure signature;
- the reacting 2D hotspot produces OH and H2O, a finite pressure/velocity response, and roundoff-scale composition closure.

The oblique wave is an exact advected entropy wave. The Gaussian velocity vortex is only a structural regression and is not presented as an analytic steady solution. The current CTU claim covers characteristic PLM normal prediction and conservative transverse correction; it does not cover multidimensional PeleC PPM tracing, physical boundaries, viscosity, heat conduction, or species diffusion.
