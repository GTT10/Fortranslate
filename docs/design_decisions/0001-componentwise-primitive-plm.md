# Decision 0001: introduce componentwise primitive PLM before PeleC characteristic PLM

## Status

Accepted as an intermediate implementation step.

## Context

PeleC's production Godunov path combines reconstruction, characteristic information, tracing, Riemann solves, and multidimensional corrections. Porting all of those pieces simultaneously would make failures difficult to localize and would provide no trustworthy intermediate parity gate.

## Decision

Implement a smaller second-order path first:

- reconstruct `rho`, velocity components, and pressure componentwise;
- provide minmod and monotonized-central limiters;
- use the existing Rusanov flux;
- retain the piecewise-constant path;
- verify linear exactness, smooth second-order convergence, shock robustness, positivity, and conservation.

## Consequences

Positive consequences:

- reconstruction infrastructure is tested independently of a new Riemann solver;
- periodic boundary handling is verified by convergence rather than visual inspection;
- the first-order path remains a diagnostic fallback;
- future characteristic PLM can be compared against a working second-order baseline.

Limitations:

- this is not a line-by-line translation of `Source/PLM.H`;
- characteristic wave coupling and PeleC-specific tracing are absent;
- matching the current Sod result does not establish full PeleC Godunov parity.

The next step is to introduce the PeleC-style approximate Riemann solver as another selectable component, then add characteristic reconstruction and tracing behind separate tests.
