# Decision 0008: generate a qualified elementary mechanism before adding stiff chemistry

## Status

Accepted for milestone `0.9.0`.

## Context

The `0.8.0` reactor verified NASA7 units, formation-energy offsets, and stage-wise adiabatic temperature recovery using a synthetic `A -> B` source. It did not test arbitrary stoichiometry, reverse rates, equilibrium constants, or production-rate parity against an external chemistry package.

Moving directly to a complete H2/O2 or hydrocarbon mechanism would simultaneously introduce mechanism parsing, third-body and falloff forms, stiffness, Jacobians, and solver coupling. A failure would be difficult to localize.

## Decision

Introduce a qualified intermediate layer:

- define a runtime elementary-reaction record with reactant/product stoichiometric vectors;
- evaluate forward Arrhenius rates in SI units;
- derive reverse rates from NASA7 equilibrium constants;
- assemble progress rates and species production rates generically;
- author a normalized JSON mechanism and generate a committed Fortran module;
- require regeneration to reproduce the committed source exactly;
- select four reversible elementary reactions from Cantera `h2o2.yaml`;
- integrate them in an adaptive explicit constant-volume reactor;
- compare both the reactor trajectory and exact-state production rates with Cantera 3.2.

## Consequences

Benefits:

- stoichiometry, units, equilibrium, source assembly, and energy coupling are independently testable;
- the generated mechanism path is reproducible rather than manually copied into Fortran;
- exact-state production-rate comparison separates kernel parity from integration error;
- Cantera parity is quantified before introducing stiff integration.

Limitations:

- the normalized JSON is not yet a Cantera YAML or CHEMKIN parser;
- third-body efficiencies, pressure falloff, and Troe/SRI are absent;
- the four-reaction system is not the complete Cantera H2/O2 mechanism;
- explicit RK4 is not suitable for a general stiff combustion mechanism;
- chemistry is not coupled to the flow solver.

The next step is third-body/falloff support, a generated Jacobian, and CVODE or equivalent stiff integration, followed by a complete small H2/O2 mechanism.
