# Decision 0006: preserve the five-equation hydro ABI behind a dynamic multispecies state

## Status

Accepted for the non-reacting multispecies milestone.

## Context

PeleC stores internal energy, temperature, and species densities in a larger state vector, while the verified PeleF Riemann and characteristic kernels currently operate on five Euler conserved variables. Expanding every existing routine to a compile-time maximum species count would obscure the tested hydro subset and create unnecessary coupling before mixture thermodynamics exists.

Species closure also must not be repaired by silently renormalizing cell averages, because doing so would hide errors in conservative flux construction.

## Decision

Use a dynamic wrapper state

```text
[rho, rho*u, rho*v, rho*w, rho*E, rho*e, T, rho*Y_1, ..., rho*Y_N]
```

with these rules:

- the first five components remain the unchanged hydro ABI;
- `rho*e` and `T` are synchronized derived fields for now;
- species count is selected at runtime from 1 to 32;
- hydro Riemann solvers return the mass flux;
- species fluxes use donor mass fractions and sum exactly to the mass flux;
- face mass fractions may be bounded and normalized;
- cell-centered species densities are never renormalized;
- 2D passive transverse corrections are applied conservatively to `rho*Y_k` using the same hydro face correction factor;
- a face-only first-order fallback is counted and exposed when a corrected passive face is invalid.

## Consequences

Benefits:

- all existing five-equation tests remain unchanged;
- species indexing and closure are independently testable;
- future mechanism sizes do not require recompiling fixed species dimensions into every hydro API;
- 1D and 2D species conservation share the same invariant;
- physical mixture thermodynamics can replace the derived placeholder without changing state positions.

Limitations:

- composition currently has no effect on pressure, sound speed, energy, or temperature;
- no molecular weights, NASA polynomials, transport, diffusion, or reactions are present;
- the 2D multispecies path is currently exercised through regression programs rather than a standalone application driver.

The next step is mechanism metadata and mixture thermodynamics, followed by temperature inversion and a non-reacting multispecies shock-tube comparison.
