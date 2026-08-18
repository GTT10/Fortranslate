# Numerical methods

## Governing equations

PeleF advances the compressible Euler system

\[
\partial_t U + \partial_x F(U) + \partial_y G(U)=0,
\qquad
U=(\rho,\rho u,\rho v,\rho w,\rho E)^T,
\]

with a constant-`gamma` ideal-gas closure. The y term is omitted by the one-dimensional driver.

## One-dimensional Godunov components

The 1D implementation retains PCM, componentwise primitive PLM, and a qualified PeleC-style characteristic PLM path. The characteristic path supports:

- `u-c`, `u`, `u+c` wave tracing;
- order-2 or PeleC five-point order-4 limited slopes;
- optional pressure/normal-velocity shock flattening;
- Rusanov or the qualified single-species PeleC-style Riemann solver.

For order 4, neighboring limited slopes enter

\[
d_{\mathrm{temp}}=
\frac{4}{3}d_{\mathrm{cen}}
-\frac{1}{6}\left(d_{i-1}+d_{i+1}\right),
\]

and the final slope is bounded by the local monotonicity constraint and multiplied by the flattening coefficient.

## Directional flux rotation

The two-dimensional solver reuses the x-direction Riemann interface for y faces. A y-normal state is rotated by exchanging x and y momentum,

\[
(\rho,\rho u,\rho v,\rho w,\rho E)
\mapsto
(\rho,\rho v,\rho u,\rho w,\rho E),
\]

then passed to the x solver. The returned momentum-flux components are exchanged again. This gives a single tested Riemann implementation per solver rather than duplicated directional code.

## Two-dimensional normal prediction

For each cell, limited primitive slopes are computed independently in x and y. Density and pressure slopes are reduced if either extrapolated face approaches its configured floor.

The x slope is projected and traced with the existing characteristic kernel. For y, primitive velocity components are rotated, traced with the same kernel using `dt/dy`, and rotated back. These predictors contain the normal half-time evolution but not the transverse flux divergence.

## CTU-style transverse correction

Let `F_hat` and `G_hat` be provisional Riemann fluxes from the normal predictor states. The left state on an x face is corrected with the transverse flux divergence of its originating cell,

\[
U^{*}_{L,i+1/2,j}
=
U^{n+1/2}_{L,i+1/2,j}
-
\frac{\Delta t}{2\Delta y}
\left(
\widehat G_{i,j+1/2}-\widehat G_{i,j-1/2}
\right).
\]

The right x-face state uses the corresponding divergence in cell `i+1`. The lower and upper states on a y face are corrected analogously with provisional x fluxes,

\[
U^{*}_{B,i,j+1/2}
=
U^{n+1/2}_{B,i,j+1/2}
-
\frac{\Delta t}{2\Delta x}
\left(
\widehat F_{i+1/2,j}-\widehat F_{i-1/2,j}
\right).
\]

Final Riemann fluxes are evaluated from the corrected states. The conservative update is

\[
U^{n+1}_{i,j}=U^n_{i,j}
-
\frac{\Delta t}{\Delta x}
\left(F_{i+1/2,j}-F_{i-1/2,j}\right)
-
\frac{\Delta t}{\Delta y}
\left(G_{i,j+1/2}-G_{i,j-1/2}\right).
\]

This is a regular-grid two-dimensional CTU-style subset. It does not yet include PeleC source-term coupling, embedded boundaries, 3D double-transverse terms, or AMR synchronization.

## Positivity scaling of transverse corrections

If a full transverse correction would produce non-positive density or pressure, PeleF seeks

\[
U(\theta)=U_{\mathrm{base}}+\theta\,\Delta U,
\qquad 0\leq\theta\le1,
\]

and accepts the largest physical `theta` found by bisection. Smooth-vortex tests accept `theta=1` everywhere.

## Two-dimensional CFL condition

The unsplit timestep uses

\[
\Delta t=
\frac{\mathrm{CFL}}
{\max_{i,j}\left[
(|u|+c)/\Delta x+(|v|+c)/\Delta y
\right]}.
\]

## Isentropic vortex

The periodic analytical test superimposes a vortex of strength `beta` on a uniform flow. With displacement `(x_bar,y_bar)` from the translated periodic vortex center,

\[
\delta u=-\frac{\beta}{2\pi}y_{\!bar}
\exp\left(\frac{1-r^2}{2}\right),
\qquad
\delta v=\frac{\beta}{2\pi}x_{\!bar}
\exp\left(\frac{1-r^2}{2}\right),
\]

\[
\delta T=-\frac{(\gamma-1)\beta^2}{8\gamma\pi^2}
\exp(1-r^2).
\]

Density and pressure follow the isentropic relation. The exact solution is the initial vortex translated by the base velocity through the periodic domain.

## Verified multidimensional results

With the PeleC-style Riemann solver, MC slopes, and transverse correction enabled:

```text
24 x 24   density L1 = 2.5862222302e-3
48 x 48   density L1 = 5.3334803639e-4
96 x 96   density L1 = 1.1010295208e-4

observed order:
24 -> 48  2.2776971
48 -> 96  2.2762241
```

At `48 x 48`, disabling the transverse correction increases density L1 error to `1.1493796865e-3`. Periodic mass, both momentum components, and total energy remain conserved to approximately `6e-14`.

## Scope limitations

The current 2D path is serial, periodic, single-species, inviscid, constant-`gamma`, and uniform-grid. General-EOS internal-energy characteristics, species/passive-scalar transport, physical wall/inflow boundaries, PPM/WENO, embedded boundaries, AMR, MPI, diffusion, chemistry, and spray remain future work.

## Passive multispecies transport

For each species, PeleF advances `rho*Y_k` conservatively. The interface flux is

```text
F_(rho Y_k) = F_rho * Y_k^upwind
```

with the final species used as a deterministic closure component so the sum of species fluxes equals the mass flux to roundoff. In 1D, mass fractions are traced with the contact-wave velocity. In 2D, species face states receive the same CTU transverse half-step correction as the hydro face states. Cell updates are accepted only when species densities remain non-negative and their sum matches total density.

## NASA7 species thermodynamics

For each species, the active seven-coefficient interval gives

\[
\frac{c_p^\circ}{R}=a_1+a_2T+a_3T^2+a_4T^3+a_5T^4,
\]

\[
\frac{h^\circ}{RT}=a_1+\frac{a_2T}{2}+\frac{a_3T^2}{3}
+\frac{a_4T^3}{4}+\frac{a_5T^4}{5}+\frac{a_6}{T},
\]

\[
\frac{s^\circ}{R}=a_1\ln T+a_2T+\frac{a_3T^2}{2}
+\frac{a_4T^3}{3}+\frac{a_5T^4}{4}+a_7.
\]

Molar values are converted to mass-specific values with the species molecular weight. Internal energy and constant-volume heat capacity follow

\[
u=h-R_kT,\qquad c_v=c_p-R_k.
\]

The returned entropy is the mass-weighted standard-state species contribution; a pressure-dependent ideal-mixing entropy term is not yet included.

## Ideal-gas mixture properties

For mass fractions `Y_k`,

\[
\frac{1}{W_{mix}}=\sum_k\frac{Y_k}{W_k},
\qquad R_{mix}=\frac{R_u}{W_{mix}}.
\]

Mass-specific caloric properties are weighted directly,

\[
c_p=\sum_kY_kc_{p,k},\quad
c_v=\sum_kY_kc_{v,k},\quad
h=\sum_kY_kh_k,\quad
u=\sum_kY_ku_k,
\]

with `gamma = cp/cv`. Pressure, density, and the frozen-composition sound speed are

\[
p=\rho R_{mix}T,
\qquad
\rho=\frac{p}{R_{mix}T},
\qquad
c=\sqrt{\gamma R_{mix}T}.
\]

These functions are currently verified independently and are not yet used by the hydro Riemann solvers.

## Internal-energy temperature inversion

The inversion solves

\[
f(T)=u(Y,T)-u_{target}=0
\]

inside the common NASA7 validity interval. The algorithm first evaluates both bracket endpoints. Each iteration attempts

\[
T_{new}=T-\frac{f(T)}{c_v(Y,T)}.
\]

A non-finite or out-of-bracket Newton candidate is replaced with the bracket midpoint. Targets outside the endpoint energy interval are rejected.

## Toy constant-volume reactor

The verification reaction is

```text
A -> B
```

with equal molecular weights and first-order rate

\[
\dot Y_A=-k(T)Y_A,\qquad
\dot Y_B=+k(T)Y_A,
\]

\[
k(T)=AT^b\exp(-T_a/T).
\]

Composition is integrated with classical RK4. In adiabatic mode, every RK stage solves

\[
u(Y^{stage},T^{stage})=u(Y^0,T^0),
\]

so the rate and heat release use a stage-consistent temperature. In isothermal mode, the analytical solution `Y_A(t)=Y_A(0) exp(-kt)` is used as a unit gate.
