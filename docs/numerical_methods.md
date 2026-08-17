# Numerical methods

## Governing equations

The current solver advances the one-dimensional Euler system in conservative form,

\[
\partial_t U + \partial_x F(U)=0,
\]

with

\[
U=(\rho,\rho u,\rho v,\rho w,\rho E)^T.
\]

The thermodynamic closure is a constant-`gamma` ideal gas.

## Finite-volume operator

For cell `i`,

\[
\frac{dU_i}{dt}=-\frac{F_{i+1/2}-F_{i-1/2}}{\Delta x}.
\]

The interface flux is Rusanov/local Lax–Friedrichs,

\[
F^*=\frac{F(U_L)+F(U_R)}{2}
-\frac{1}{2}\max(|u_L|+c_L,|u_R|+c_R)(U_R-U_L).
\]

## Reconstruction choices

### Piecewise constant

`reconstruction = "pcm"` sends adjacent cell averages directly to the Riemann solver. This path remains permanently available as a robustness and regression baseline.

### Piecewise linear

`reconstruction = "plm"` converts cell averages to primitive variables

\[
q=(\rho,u,v,w,p)^T,
\]

computes a limited componentwise slope `s_i`, and forms

\[
q_{i+1/2}^{L}=q_i+\frac{1}{2}s_i,
\qquad
q_{i-1/2}^{R}=q_i-\frac{1}{2}s_i.
\]

Face primitive states are converted back to conserved variables before flux evaluation.

## Slope limiters

For backward and forward differences `a` and `b`, minmod is

\[
\operatorname{minmod}(a,b)=
\begin{cases}
\operatorname{sign}(a)\min(|a|,|b|), & ab>0,\\
0, & \text{otherwise}.
\end{cases}
\]

The monotonized-central slope is

\[
s_i=\operatorname{minmod}\left(
\frac{a+b}{2},2a,2b
\right).
\]

Available input values are `minmod` and `mc`.

## Physical-state protection

Density and pressure slopes share a scalar reduction factor if either extrapolated face would reach its configured floor. If conversion of an extrapolated primitive state still fails, that face reverts to the corresponding cell-centered conserved state. The cell update itself is rejected if SSPRK2 produces a non-physical interior state.

## Boundary conditions

- `outflow`: nearest interior state copied into the ghost cell; slopes at boundary-adjacent cells are suppressed.
- `periodic`: states and slopes wrap between the first and last cells.

## Time integration

The semi-discrete operator is advanced by SSPRK2:

\[
U^{(1)}=U^n+\Delta t L(U^n),
\]

\[
U^{n+1}=\frac{1}{2}U^n+
\frac{1}{2}\left(U^{(1)}+\Delta t L(U^{(1)})\right).
\]

The timestep is limited by

\[
\Delta t=\mathrm{CFL}\min_i\frac{\Delta x}{|u_i|+c_i}.
\]

## Current accuracy evidence

The periodic entropy-wave regression produces density L1 errors of `4.8297e-4`, `1.1659e-4`, and `2.7383e-5` on 40, 80, and 160 cells. The corresponding observed orders are `2.0505` and `2.0901`.

For the 400-cell Sod case at `t = 0.2`, PLM/MC reduces the density L1 error from approximately `1.195e-2` to `1.891e-3` and pressure L1 from approximately `9.244e-3` to `1.198e-3` relative to the retained PCM baseline.

## Scope limitation

This PLM implementation is componentwise in primitive variables. PeleC's characteristic projection, characteristic tracing, flattening details, and multidimensional corrections are not yet implemented and must be verified independently before algorithmic parity is claimed.
