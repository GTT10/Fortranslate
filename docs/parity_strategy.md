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

## Chemistry-kernel gates

The shared reaction layer verifies:

- SI Arrhenius evaluation;
- valid and invalid stoichiometric records;
- reversible equilibrium constants from NASA7 Gibbs functions;
- elementary, third-body, and falloff progress rates;
- collider efficiencies and Troe broadening;
- production-rate and mass-fraction-source assembly;
- instantaneous mass and H/O atom conservation;
- analytic Jacobians against independent centered finite differences.

Both generated H2/O2 modules must regenerate byte-for-byte from their normalized JSON inputs.

## Reactor structural gates

The application-level checkers read only emitted CSV files and independently enforce:

- strictly increasing output times and the requested final time;
- constant density;
- finite, non-negative mass fractions and closure;
- fixed specific internal energy;
- H and O atom inventories;
- unchanged inert species;
- instantaneous mass and atom conservation of production rates;
- non-trivial temperature and composition evolution.

The implicit unit gates also require Newton convergence, physical line-search states, a positive accepted step, and agreement of the reduced constant-energy Jacobian with finite differences.

## Live Cantera parity

CI installs Cantera 3.2.0 and maintains two independent gates.

### Elementary subset

The seven-species/four-reaction companion YAML is compared at tight tolerances. This remains the fast high-precision rate and explicit-integrator regression.

### Full pressure-dependent mechanism

Cantera loads its pinned `h2o2.yaml` `ohmech` phase. The PeleF and Cantera reactors begin from the same temperature, pressure, and ten-species composition and are sampled at the same 101 output times. Temperature, pressure, and every species mass fraction are compared.

At every PeleF output row, a separate Cantera phase is reset to the exact PeleF `(T,rho,Y)` state and its `net_production_rates` are compared with the Fortran kernel. This separates kinetic-kernel error from integration-history error.

Current full-mechanism trajectory differences are approximately:

```text
maximum temperature difference      4.59e-3 K
maximum pressure difference         3.72e-1 Pa
maximum species relative difference < 1.0e-5
final-temperature difference        7.18e-5 K
```

The full trajectory relative tolerances are `5e-5` for temperature/pressure and `2e-4` for species. Production rates use `2e-8` relative tolerance and a `1e-10 kmol/(m^3 s)` absolute floor for nearly cancelled net sources. Thresholds may change only with an explained units, mechanism-data, or numerical-method change.

## Reference-data policy

Every external comparison records:

- upstream repository and commit SHA;
- source mechanism and phase;
- units and SI conversions;
- species ordering and molecular weights;
- initial state and reactor model;
- solver tolerances;
- comparison variables and sample times.

Pinned numerical signatures may be updated only with an explained method or data change. Conservation limits are never relaxed merely to accept a regression.

## Scope of the evidence

The current evidence establishes zero-dimensional parity for the complete small Cantera H2/O2 mechanism and the implemented elementary/third-body/Troe forms. It does not establish general PelePhysics parser parity, SRI/chemically activated reactions, production-scale sparse stiff integration, hydrocarbon chemistry, or chemistry-coupled CFD.
