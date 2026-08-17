# Numerical methods

## Governing equations

The current solver advances the one-dimensional Euler system in conservative form,

\[
\partial_t U + \partial_x F(U)=0,
\qquad
U=(\rho,\rho u,\rho v,\rho w,\rho E)^T,
\]

with a constant-`gamma` ideal-gas closure.

For cell `i`, the finite-volume operator is

\[
\frac{dU_i}{dt}=-\frac{F_{i+1/2}-F_{i-1/2}}{\Delta x}.
\]

## Reconstruction

`reconstruction = "pcm"` sends adjacent cell averages directly to the selected Riemann solver.

`reconstruction = "plm"` converts cell averages to primitive variables

\[
q=(\rho,u,v,w,p)^T,
\]

computes a componentwise limited slope `s_i`, and forms

\[
q_{i+1/2}^{L}=q_i+\frac{1}{2}s_i,
\qquad
q_{i-1/2}^{R}=q_i-\frac{1}{2}s_i.
\]

Available limiters are minmod and monotonized central (`mc`). Density and pressure slopes share a reduction factor if an extrapolated face approaches a configured physical floor. If primitive-to-conserved conversion still fails, the affected face falls back to the corresponding cell-centered state.

## Riemann solver choices

### Rusanov baseline

The robust baseline is local Lax-Friedrichs/Rusanov,

\[
F^*=\frac{F(U_L)+F(U_R)}{2}
-\frac{1}{2}\max(|u_L|+c_L,|u_R|+c_R)(U_R-U_L).
\]

It remains selectable for differential diagnosis and future difficult-state tests.

### PeleC-style approximate solver

For the current ideal-gas, single-species subset, the acoustic impedances are

\[
Z_L=\rho_Lc_L,\qquad Z_R=\rho_Rc_R.
\]

The initial star estimates follow the corresponding PeleC `Source/Riemann.H` structure,

\[
p^*=\frac{Z_Rp_L+Z_Lp_R+Z_LZ_R(u_L-u_R)}{Z_L+Z_R},
\]

\[
u^*=\frac{Z_Lu_L+Z_Ru_R+p_L-p_R}{Z_L+Z_R}.
\]

The origin state is selected by the sign of `u*`, with stationary interfaces averaged. Star density is estimated by

\[
\rho^*=\rho_o+\frac{p^*-p_o}{c_o^2}.
\]

PeleC-style inward and outward wave speeds then interpolate between the origin and star states, with shock replacement when `p* >= p_o`. The resulting interface state supplies mass, momentum, and total-energy fluxes.

This implementation intentionally excludes the production solver's general PelePhysics EOS, species arrays, rotating-frame terms, and boundary scaling. It is therefore described as a verified reduction, not full `Source/Riemann.H` parity.

## Time integration

The semi-discrete operator is advanced by SSPRK2,

\[
U^{(1)}=U^n+\Delta t L(U^n),
\]

\[
U^{n+1}=\frac{1}{2}U^n+
\frac{1}{2}\left(U^{(1)}+\Delta t L(U^{(1)})\right),
\]

with

\[
\Delta t=\mathrm{CFL}\min_i\frac{\Delta x}{|u_i|+c_i}.
\]

## Current accuracy evidence

For a periodic entropy wave using PLM/MC and the PeleC-style solver, density L1 errors on 40, 80, and 160 cells are approximately `4.6967e-4`, `1.1084e-4`, and `2.5946e-5`; observed orders are `2.0831` and `2.0949`.

For the 400-cell Sod case at `t = 0.2`:

| Method | Density L1 | Pressure L1 |
|---|---:|---:|
| PCM + Rusanov | `1.1952e-2` | `9.2441e-3` |
| PLM/MC + Rusanov | `1.8907e-3` | `1.1981e-3` |
| PLM/MC + PeleC-style | `1.3678e-3` | `8.0552e-4` |

The Shu-Osher case reaches `t = 1.8` on 800 cells with positive density and pressure, roundoff-scale integral-balance errors, and 21 detected extrema in the interaction window.

## Remaining Godunov work

Characteristic-variable projection, characteristic tracing, flattening, multidimensional transverse corrections, PPM, and WENO remain unimplemented. Those are separate parity gates rather than implied by the Riemann-solver result.
