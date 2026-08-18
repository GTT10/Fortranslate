# Decision 0007: verify NASA7 thermodynamics before hydro coupling

## Status

Accepted for milestone `0.8.0`.

## Context

The passive multispecies solver transports `rho*Y_k`, but its hydro closure remains constant-`gamma`. Replacing that closure while simultaneously adding molecular weights, NASA polynomials, temperature inversion, and chemistry would make any regression difficult to localize.

## Decision

Add an independent thermodynamics and zero-dimensional reactor layer first:

- store molecular weight, validity bounds, and low/high NASA7 coefficients per species;
- evaluate species caloric properties in SI mass units;
- form ideal-gas mixture properties from mass fractions;
- solve `u(Y,T)=u_target` with a bracketed Newton/bisection method;
- verify H2, O2, H2O, and N2 data from the public Cantera mechanism files;
- test source integration with a two-species equal-molecular-weight toy isomerization;
- keep all existing hydro Riemann and CFL paths unchanged.

## Consequences

Benefits:

- coefficient evaluation, units, mixture algebra, and temperature inversion have isolated tests;
- formation-energy offsets are exercised by the adiabatic reactor;
- the reactor can verify stage-consistent thermodynamic coupling without a large mechanism;
- existing constant-`gamma` hydro signatures remain valid.

Limitations:

- the toy reaction is not physical combustion chemistry;
- no mechanism parser, reverse rate, Jacobian, or stiff solver exists;
- standard-state entropy is available, but ideal-mixing entropy is not yet assembled;
- NASA7 thermodynamics are not yet wired into hydro pressure, sound speed, or characteristic tracing.

The next step is an arbitrary elementary-reaction representation and generated production-rate kernel, followed by direct Cantera comparison on a real small mechanism.
