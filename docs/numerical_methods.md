# Numerical methods

## Governing equations

PeleF advances the compressible Euler system

\[
\partial_t U + \partial_x F(U) + \partial_y G(U)=0,
\qquad
U=(\rho,\rho u,\rho v,\rho w,\rho E)^T,
\]

The established `pelef` and `pelef2d` paths use a constant-`gamma` ideal-gas closure. The separate `pelef_reactive_1d` path appends species densities and closes pressure, temperature, and sound speed with NASA7 ideal-gas-mixture thermodynamics.

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

## Generated mechanism subset

`tools/generate_elementary_mechanism.py` reads a normalized JSON mechanism and writes a Fortran module containing fixed species indices, reaction stoichiometry, SI Arrhenius parameters, and a production-rate wrapper. CI regenerates the module and requires exact agreement with the committed source.

The current H2/O2 subset contains four reversible elementary reactions:

```text
O + H2  <=> H + OH
H + O2  <=> O + OH
OH + H2 <=> H + H2O
2 OH    <=> O + H2O
```

H2, H, O, O2, OH, H2O, and N2 are present. N2 is inert because no third-body reactions are included yet.

## Adaptive constant-volume reactor

The H2/O2 reactor uses explicit RK4 with step doubling:

1. one full RK4 step of size `dt`;
2. two RK4 half steps;
3. a scaled norm of the state difference estimates local error;
4. the two-half-step result is accepted when the error is within tolerance;
5. the next step is expanded or contracted within configured bounds.

At every RK stage, temperature is recovered from

\[
u(Y^{\mathrm{stage}},T^{\mathrm{stage}})=u(Y^0,T^0),
\]

so the rate evaluation and adiabatic energy constraint are stage consistent. Trial states are rejected when composition loses positivity or closure.

This explicit integrator is suitable only for the current verification subset. A stiff integrator and Jacobian are required before using a complete combustion mechanism.

## Verified H2/O2 parity

The live Cantera gate separates two comparisons:

- trajectory parity: temperature, pressure, and all species at each output time;
- kinetic-kernel parity: Cantera production rates evaluated at the exact PeleF `(T,rho,Y)` state.

For the current case, the maximum absolute differences are:

```text
temperature              1.61e-6 K
pressure                 1.36e-4 Pa
species mass fraction    1.70e-11
production rate          3.55e-12 kmol/(m^3 s)
final temperature        3.69e-9 K
```

The production-rate maximum is an almost cancelled OH net source near `2.5e-8 kmol/(m^3 s)`, so parity uses both a relative tolerance and a `5e-12 kmol/(m^3 s)` absolute floor.

## General-EOS reactive conserved state

The reactive one-dimensional driver advances

\[
U=\left(\rho,\rho u,\rho v,\rho w,\rho E,
\rho Y_1,\ldots,\rho Y_N\right)^T.
\]

For each cell, the specific internal energy is recovered from

\[
e=\frac{\rho E}{\rho}
-\frac{1}{2}\left(u^2+v^2+w^2\right).
\]

The nonlinear equation

\[
u(Y,T)=e
\]

is then solved in the shared NASA7 validity interval. Pressure and frozen-composition sound speed follow from

\[
p=\rho R_{\mathrm{mix}}(Y)T,
\qquad
c=\sqrt{\gamma(Y,T)\frac{p}{\rho}}.
\]

The inverse conversion reconstructs total energy with the full NASA7 internal energy, including formation-energy offsets.

## Frozen-composition characteristic PLM

The current reactive high-order path reconstructs

\[
q=(\rho,u,v,w,p,Y_1,\ldots,Y_N)^T.
\]

For the hydrodynamic prefix, differences are projected onto five frozen-composition waves with speeds

```text
u-c, u, u, u, u+c.
```

The acoustic and contact strengths use the local mixture sound speed:

\[
\alpha_- = \frac{\rho}{2c}
\left(\frac{\Delta p}{\rho c}-\Delta u\right),
\qquad
\alpha_+ = \frac{\rho}{2c}
\left(\frac{\Delta p}{\rho c}+\Delta u\right),
\]

\[
\alpha_0=\Delta\rho-\frac{\Delta p}{c^2}.
\]

Each characteristic difference is limited with minmod or MC. Species mass fractions are limited componentwise and traced at the contact-wave speed `u`. A MUSCL-Hancock prediction forms time-centered face states,

\[
q_{L/R}^{n+1/2}
=q_i\mp\frac{1}{2}\Delta q_i
-\frac{\Delta t}{2\Delta x}A_i\Delta q_i.
\]

Density, pressure, and mass-fraction slopes are scaled before tracing when an unlimited face would leave the physical interval. Face mass fractions are normalized, and invalid faces fall back to the cell-centered state.

This is a qualified ideal-gas-mixture, frozen-composition basis. It is not yet the complete general-EOS characteristic treatment used through PelePhysics.

## Reactive Rusanov flux and timestep

The physical flux is

\[
F(U)=\left(
\rho u,
\rho u^2+p,
\rho uv,
\rho uw,
(\rho E+p)u,
\rho uY_k
\right)^T.
\]

The interface flux is

\[
\widehat F=rac{F_L+F_R}{2}
-\frac{a_{\max}}{2}(U_R-U_L),
\qquad
a_{\max}=\max(|u_L|+c_L,|u_R|+c_R).
\]

The last species flux is assigned from the total mass flux minus the preceding species fluxes, enforcing

\[
\sum_k F_{\rho Y_k}=F_\rho.
\]

The timestep uses the maximum composition-dependent signal speed,

\[
\Delta t=\mathrm{CFL}
\frac{\Delta x}{\max_i(|u_i|+c_i)}.
\]

## Strang-split reaction-flow coupling

A complete step is

\[
U^n
\xrightarrow{R(\Delta t/2)}U^{(1)}
\xrightarrow{H(\Delta t)}U^{(2)}
\xrightarrow{R(\Delta t/2)}U^{n+1}.
\]

During `R`, each cell retains `rho`, all momenta, and `rhoE`. The constant-volume reactor changes `Y`; temperature is recovered from the unchanged internal energy after every accepted chemistry substep. During `H`, all conserved components, including `rho*Y_k`, are updated by the finite-volume flux divergence.

A homogeneous periodic field is required to reduce to the independent zero-dimensional reactor. A nonuniform temperature hotspot tests the coupled generation and propagation of pressure and velocity disturbances.

## Reactive-flow numerical evidence

The periodic pure-N2 entropy wave on 40, 80, and 160 cells gives density L1 errors

```text
1.51594309e-4
3.43297270e-5
7.01334896e-6
```

and observed orders `2.142685` and `2.291283`.

For the 30 microsecond H2/O2 hotspot regression, the 64-cell solution remains positive and produces

```text
pressure span       1.0382765e3 Pa
maximum |u|         2.2117421e1 m/s
temperature range   1328.97 to 1545.46 K
```

A separate 128-cell characteristic-PLM solution is used as a restricted reference at 8 microseconds. The combined normalized L1 errors are

```text
                 32 cells       64 cells
characteristic PLM  6.2286e-2      1.4690e-2
PCM                 2.0249e-1      1.5436e-1
```

so refinement reduces the PLM error, and PLM materially outperforms the first-order baseline at both resolutions.

## Scope limitations

The code remains serial and uniform-grid. Reactive CFD currently uses only the seven-species, four-reaction elementary subset, Rusanov fluxes, and a frozen-composition characteristic basis. It lacks third-body/falloff chemistry, a stiff cell integrator, species diffusion, viscosity, thermal conduction, a complete mechanism, general-EOS PeleC Riemann parity, multidimensional reactive flow, AMR, MPI, and accelerators.
