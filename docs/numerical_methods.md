# Numerical methods

## Governing equations

PeleF advances the compressible Euler system

\[
\partial_t U + \partial_x F(U) + \partial_y G(U)=0,
\qquad
U=(\rho,\rho u,\rho v,\rho w,\rho E)^T,
\]

with a constant-`gamma` ideal-gas closure. The y term is omitted by the one-dimensional driver. The current NASA7 and chemistry layers are not yet used by these fluxes.

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

## Two-dimensional CTU-style update

The 2D solver rotates y-normal states into the x-normal Riemann interface, computes independent normal predictions, evaluates provisional fluxes, applies transverse half-step flux differences, recomputes final face fluxes, and performs one unsplit conservative update.

For example, the left state on an x face is corrected by

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

If a full correction would produce non-positive density or pressure, the correction is scaled by the largest physical `0 <= theta <= 1` found by bisection.

The 2D timestep is

\[
\Delta t=
\frac{\mathrm{CFL}}
{\max_{i,j}\left[(|u|+c)/\Delta x+(|v|+c)/\Delta y\right]}.
\]

## Passive multispecies transport

For each passive species, PeleF advances `rho*Y_k` conservatively. The interface flux is

\[
F_{\rho Y_k}=F_\rho Y_k^{\mathrm{upwind}},
\]

with a deterministic closure component so that the species flux sum equals the mass flux to roundoff. In 1D, mass fractions are traced with the contact-wave velocity. In 2D, species face states receive the same CTU transverse correction as the hydro face states.

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

## Ideal-gas mixture properties

For mass fractions `Y_k`,

\[
\frac{1}{W_{\mathrm{mix}}}=\sum_k\frac{Y_k}{W_k},
\qquad R_{\mathrm{mix}}=\frac{R_u}{W_{\mathrm{mix}}}.
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
p=\rho R_{\mathrm{mix}}T,
\qquad
\rho=\frac{p}{R_{\mathrm{mix}}T},
\qquad
c=\sqrt{\gamma R_{\mathrm{mix}}T}.
\]

The temperature inversion solves `u(Y,T)=u_target` inside the common NASA7 validity interval using Newton updates with mixture `cv` and a bisection fallback.

## Elementary reaction representation

Each reaction stores reactant and product stoichiometric vectors,

\[
\nu'_{k,r},\qquad \nu''_{k,r},
\]

an Arrhenius triplet `(A,b,E_a)`, and a reversible flag. The forward rate constant is

\[
k_{f,r}=A_rT^{b_r}\exp\left(-\frac{E_{a,r}}{R_uT}\right).
\]

Molar concentrations in `kmol/m^3` are recovered from density and mass fractions:

\[
C_k=\frac{\rho Y_k}{W_k}.
\]

The forward progress rate is

\[
q_{f,r}=k_{f,r}\prod_k C_k^{\nu'_{k,r}}.
\]

## Reverse rates and equilibrium constants

NASA7 standard-state Gibbs functions give

\[
\frac{\Delta g_r^\circ}{R_uT}
=
\sum_k(\nu''_{k,r}-\nu'_{k,r})
\left(\frac{h_k^\circ}{R_uT}-\frac{s_k^\circ}{R_u}\right).
\]

With `Delta nu_r = sum_k(nu''-nu')`, the concentration equilibrium constant is

\[
K_{c,r}
=
\exp\left(-\frac{\Delta g_r^\circ}{R_uT}\right)
\left(\frac{p^\circ}{R_uT}\right)^{\Delta\nu_r}.
\]

The reverse constant and progress rate are

\[
k_{r,r}=\frac{k_{f,r}}{K_{c,r}},
\qquad
q_{r,r}=k_{r,r}\prod_k C_k^{\nu''_{k,r}}.
\]

The net progress and species production rates are

\[
q_r=q_{f,r}-q_{r,r},
\qquad
\dot\omega_k=\sum_r(\nu''_{k,r}-\nu'_{k,r})q_r,
\]

and the constant-volume composition equation is

\[
\frac{dY_k}{dt}=\frac{W_k\dot\omega_k}{\rho}.
\]

## Generated mechanism families

`tools/generate_elementary_mechanism.py` reads normalized JSON and writes a Fortran module containing fixed species indices, stoichiometry, reaction kinds, SI Arrhenius parameters, third-body efficiencies, and Troe parameters. CI regenerates both committed modules and requires byte-for-byte agreement.

The fast subset contains four reversible elementary reactions over H2, H, O, O2, OH, H2O, and N2. The full mechanism contains H2, H, O, O2, OH, H2O, HO2, H2O2, Ar, and N2 with all 29 reactions from the pinned Cantera `h2o2.yaml`, including duplicate, third-body, and Troe falloff reactions.

## Third-body and falloff rates

For a third-body reaction, the effective collider concentration is

\[
[M]_{\mathrm{eff}}=\sum_k\alpha_k C_k.
\]

For falloff reactions,

\[
P_r=\frac{k_0[M]_{\mathrm{eff}}}{k_\infty},
\qquad
k_{\mathrm{L}}=k_\infty\frac{P_r}{1+P_r}.
\]

When Troe broadening is enabled,

\[
F_{\mathrm{cent}}=(1-a)e^{-T/T_3}+ae^{-T/T_1}+e^{-T_2/T},
\]

and the standard logarithmic Troe formula produces `F`. The effective forward constant is `k_L F`. The same effective forward constant is divided by the NASA7 concentration equilibrium constant for the reverse reaction.

## Chemistry Jacobians

The production-rate kernel differentiates concentration products and pressure-dependent rate constants analytically with respect to species concentrations. Conversion to a fixed-temperature mass-fraction Jacobian uses

\[
\frac{\partial \dot Y_i}{\partial Y_j}
=\frac{W_i}{W_j}\frac{\partial \dot\omega_i}{\partial C_j}.
\]

The constant-energy reactor eliminates the final mass fraction through closure. Its reduced Jacobian includes the temperature response

\[
\frac{\partial T}{\partial Y_j}
=-\frac{u_j-u_N}{c_v},
\]

combined with a centered temperature derivative of the reaction source. The generated full-mechanism module supplies a size-checked wrapper around this shared Jacobian assembly.

## Implicit constant-volume reactor

The full mechanism solves backward Euler steps:

\[
R(Y^{n+1})=Y^{n+1}-Y^n-\Delta t\,\dot Y(Y^{n+1},T^{n+1})=0,
\]

subject to

\[
u(Y^{n+1},T^{n+1})=u_0.
\]

Newton corrections use a dense pivoted linear solve and a backtracking line search. Trial states outside the composition simplex, outside the NASA7 temperature interval, or without decreasing residual are rejected.

Adaptive control compares one full backward Euler step with two half steps. Since backward Euler is first order, the accepted candidate is

\[
Y_{\mathrm{R}}=2Y_{h/2,h/2}-Y_h.
\]

This Richardson extrapolation cancels the leading local error. When the extrapolated state is not physical, the verified two-half-step state is retained. Temperature is then recovered again from the fixed internal energy.

The seven-species subset continues to use adaptive explicit RK4. Retaining both paths separates the pressure-dependent/stiff implementation from the earlier elementary-kinetics gate.

## Verified full H2/O2 parity

The 1000 K, 1 atm, adiabatic constant-volume case is emitted at 101 fixed times through 2 ms. Against the pinned Cantera 3.2 reference, the current maximum trajectory differences are approximately:

```text
temperature              4.59e-3 K
pressure                 3.72e-1 Pa
species relative error   < 1.0e-5
final temperature        7.18e-5 K
```

The live CI gate separately resets Cantera to every exact PeleF `(T,rho,Y)` state and compares molar production rates. Near-cancelled net rates use an absolute floor of `1e-10 kmol/(m^3 s)` in addition to the relative tolerance.

## Scope limitations

The code remains serial and uniform-grid. The full chemistry path still lacks SRI and chemically activated forms, direct Cantera/CHEMKIN parsing, sparse Jacobians, CVODE/ARKODE, hydrocarbon mechanisms, molecular diffusion, and chemistry coupling to the flow solver. The dense in-tree implicit solver is a verification implementation for small mechanisms.
