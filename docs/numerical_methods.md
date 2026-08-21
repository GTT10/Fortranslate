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

## Reactive approximate Riemann fluxes and timestep

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

The robustness baseline is the local Lax--Friedrichs/Rusanov flux,

\[
\widehat F_{\mathrm{Rusanov}}=
\frac{F_L+F_R}{2}
-\frac{a_{\max}}{2}(U_R-U_L),
\qquad
a_{\max}=\max(|u_L|+c_L,|u_R|+c_R).
\]

The selectable contact-resolving path is HLLC. Conservative Davis wave bounds
are evaluated from the thermodynamic sound speed on each side,

\[
S_L=\min(u_L-c_L,u_R-c_R),
\qquad
S_R=\max(u_L+c_L,u_R+c_R).
\]

The contact speed is

\[
S_M=
\frac{
 p_R-p_L
 +\rho_Lu_L(S_L-u_L)
 -\rho_Ru_R(S_R-u_R)
}{
 \rho_L(S_L-u_L)-\rho_R(S_R-u_R)
}.
\]

For side \(K\in\{L,R\}\), the HLLC star density and energy density are

\[
\rho_K^*=\rho_K\frac{S_K-u_K}{S_K-S_M},
\]

\[
(\rho E)_K^*=
\frac{
 (S_K-u_K)(\rho E)_K-p_Ku_K+p^*S_M
}{S_K-S_M}.
\]

Tangential velocities and mass fractions remain frozen across the corresponding
outer wave. Every candidate star state is checked by the NASA7
conserved-to-primitive inversion before its flux is accepted. No silent
Rusanov fallback is used for an invalid selected HLLC state.

For both solvers, the last species flux is assigned from the total mass flux
minus the preceding species fluxes, enforcing

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

## General-EOS HLLC numerical evidence

A periodic H2/N2 composition wave couples changing molecular weight to density
while retaining uniform pressure and temperature. HLLC with characteristic PLM
gives H2 mass-fraction L1 errors

```text
40 cells     2.73882775e-5
80 cells     6.85214336e-6
160 cells    1.65453783e-6
```

with observed orders `1.998931` and `2.050127`. The corresponding 160-cell
relative density L1 error is `1.55234844e-5`.

A discontinuous moving material contact is also tested with identical pressure
and velocity but different H2/N2 composition. At 200 cells, the H2
mass-fraction L1 errors are

```text
HLLC       1.14269289e-4
Rusanov    1.87834655e-4
```

so HLLC reduces the contact error by a factor of `1.6438`. On the smooth
composition wave, Rusanov happens to have a smaller 80-cell H2 error; the tests
do not assert that HLLC is uniformly superior. The HLLC claim is limited to
contact resolution and its independently verified second-order smooth-wave
behavior.

The HLLC hotspot case remains positive and conservative while producing a
pressure span of approximately `1.3780e3 Pa` and a maximum velocity magnitude
of approximately `2.0841e1 m/s`.

## Monotone reactive PPM

The optional reactive PPM path reconstructs the general-EOS primitive vector

\[
q=(\rho,u,v,w,p,Y_1,\ldots,Y_N)^T.
\]

For a face between cells \(i\) and \(i+1\), the initial fourth-order candidate
is

\[
q_{i+1/2}=\frac{7}{12}(q_i+q_{i+1})
-\frac{1}{12}(q_{i-1}+q_{i+2}).
\]

Each component is first restricted to the range spanned by the two adjacent
cells. The left and right edge values in a cell are then constrained so that a
cell-center extremum becomes constant and excessive parabolic curvature is
removed. Density and pressure retain their physical floors; species edges are
bounded and normalized before conversion back to conserved variables.

Unlike the time-traced characteristic PLM path, this PPM reconstruction is used
as a semidiscrete spatial operator and advanced with the three-stage SSPRK3
scheme. Rusanov and HLLC remain independent flux choices. This is a monotone
primitive PPM subset.  It remains available as an independent comparison path
after addition of the separate time-traced `characteristic_ppm` algorithm.

## Time-traced reactive characteristic PPM

`reconstruction = "characteristic_ppm"` is a separate normal-predictor path;
it does not replace the semidiscrete componentwise `ppm` baseline.  For every
primitive component, the five-cell stencil is reconstructed with the van-Leer
limited edge formula used by PeleC `PPM.H`.  A cell parabola with center value
`q_c` and edge values `q_-`, `q_+` is represented through

\[
q_6 = 6q_c-3(q_-+q_+).
\]

For a characteristic speed \(\lambda\),
\(\sigma=|\lambda|\Delta t/\Delta x\).  The average swept toward the right
edge is

\[
I_+(\lambda)=
\begin{cases}
q_+, & \lambda\le 0,\\
q_+-\frac{\sigma}{2}
\left[q_+-q_--\left(1-\frac{2\sigma}{3}\right)q_6\right],
& \lambda>0,
\end{cases}
\]

and the corresponding average swept toward the left edge is

\[
I_-(\lambda)=
\begin{cases}
q_-+\frac{\sigma}{2}
\left[q_+-q_-+\left(1-\frac{2\sigma}{3}\right)q_6\right],
& \lambda\le 0,\\
q_-, & \lambda>0.
\end{cases}
\]

These integrals are evaluated for `u-c`, `u`, and `u+c`.  Species mass
fractions and transverse velocities use the middle-wave integral.  Density,
normal velocity, and pressure use the same frozen-composition acoustic/contact
projection as the characteristic PLM path, with the fastest wave moving toward
the face selected as the reference state.  Internal energy is not independently
traced in this qualified subset; the final face energy and temperature are
recovered from density, pressure, and normalized composition through the NASA7
EOS.

The reconstruction is already time-centered, so it advances with one
conservative Godunov update rather than the SSPRK3 loop used by the semidiscrete
componentwise PPM path.  A profile is rejected when the local characteristic
Courant number exceeds one.

### Reactive shock flattening

The optional `ppm_shock_flattening` detector is the one-dimensional regular-cell
formula from PeleC `Godunov.H`.  It combines a pressure-jump ratio, compression
of the normal velocity, and a shifted neighboring detector.  The resulting
coefficient

\[
f_i=1-\max(\chi_i z_i,\chi_{i+s}z_{i+s}),\qquad 0\le f_i\le1,
\]

blends each reconstructed edge toward its cell-center value before the
parabolic monotonicity constraint.  Smooth regions and expansions retain
`f=1`; the strong-compression unit stencil reaches `f=0`.

### Bounded contact steepening

The optional `ppm_contact_steepening` detector follows the
Colella--Woodward-style density/pressure criteria: a density jump must dominate
the pressure jump, neighboring density curvatures must change sign, and the
relative density change must exceed one percent.  Its canonical coefficient is

\[
\eta=\max\left[0,\min\left(1,
20(\widetilde\eta-0.05)\right)\right],
\qquad
\widetilde\eta=-\frac{\Delta^2\rho_+-\Delta^2\rho_-}
{6\Delta\rho}.
\]

Density and mass-fraction edges are blended toward neighboring MC face values
and clipped to the adjacent-cell range.  The current general-EOS HLLC subset
caps the applied strength at `0.5`; full-strength simultaneous density and
composition steepening can over-compress the material interface before a
complete PeleC/PelePhysics characteristic system is available.

### Characteristic-PPM numerical evidence

The smooth entropy-wave density L1 errors on 32, 64, and 128 cells are

```text
2.49716040e-4
4.88601959e-5
9.87993492e-6
```

with observed orders `2.353557` and `2.306086`.  The corresponding H2
mass-fraction errors for the smooth composition wave are

```text
3.31512363e-5
7.50581812e-6
1.45326368e-6
```

with observed orders `2.142981` and `2.368713`.

For the 200-cell moving material contact,

```text
componentwise PPM + HLLC             8.75871389e-5
characteristic PPM + HLLC            7.35878653e-5
characteristic PPM + bounded steepening 2.50998077e-5
```

in H2 mass-fraction L1 error.  The periodic pressure-ratio-three shock test
remains positive and conservative with or without flattening, produces no
pressure overshoot outside the initial extrema, and records a nonzero
flattened/unflattened state difference of `8.08412674e1` in its normalized
state-sum signature.

On the smooth reacting hotspot, the 32/64-cell characteristic-PPM errors against
the restricted 128-cell characteristic-PLM reference are `6.00777009e-2` and
`3.07014574e-2`.  This is convergent and far below PCM, but it is not lower than
the existing characteristic-PLM or semidiscrete componentwise-PPM error for
that particular Gaussian case.  The regression records this limitation rather
than treating the new path as universally superior.

## General-EOS reactive two-dimensional CTU

For each cell, the 2D solver recovers the primitive vector

\[
q=(\rho,u,v,w,p,Y_1,\ldots,Y_N)^T
\]

and the frozen mixture sound speed from the conserved state and synchronized
temperature guess. The normal reconstruction can be PCM, characteristic PLM, or
characteristic PPM. The y-normal calculation swaps `u` and `v`, applies the same
frozen-composition acoustic/contact basis as the x direction, and rotates the
result back.

For `characteristic_ppm`, each direction independently builds the same
five-point parabolic edge values used by the qualified 1D path. The profile is
integrated over the regions swept by `u-c`, `u`, and `u+c`; density, normal
velocity, and pressure are projected over the frozen mixture characteristic
basis, while species and transverse velocities travel on the middle wave. The
one-dimensional PeleC shock-flattening detector and the bounded contact
steepener are optional in both directions.

Provisional normal fluxes are computed at every x and y face. Each time-centered
face state then receives the transverse half-step correction

\[
U^*_{i+1/2,j}=U^{n+1/2}_{i+1/2,j}
-\frac{\Delta t}{2\Delta y}
\left(G_{i,j+1/2}-G_{i,j-1/2}\right),
\]

with the analogous x correction on y faces. The correction includes density,
all momentum components, total energy, and every species density. If the full
correction is not EOS-admissible, a scalar \(\theta\in[0,1]\) is found by
bisection and applied to the complete correction vector. The final face state
must recover positive temperature and pressure and satisfy species closure
before the final HLLC or Rusanov solve.

The cell update is unsplit:

\[
U^{n+1}_{i,j}=U^n_{i,j}
-\frac{\Delta t}{\Delta x}(F_{i+1/2,j}-F_{i-1/2,j})
-\frac{\Delta t}{\Delta y}(G_{i,j+1/2}-G_{i,j-1/2}).
\]

The timestep uses

\[
\Delta t=\mathrm{CFL}\left[\max_{i,j}\left(
\frac{|u|+c}{\Delta x}+\frac{|v|+c}{\Delta y}
\right)\right]^{-1}.
\]

Chemistry is applied cell by cell in Strang order around this hydro step. During
each reaction half-step, density, all momenta, and total energy remain fixed,
and temperature is recovered from the unchanged specific internal energy after
composition changes.

This is a normal-predictor-plus-conservative-CTU subset. It does not yet
reproduce PeleC's complete multidimensional PPM transverse/corner tracing.

### Two-dimensional numerical evidence

For the characteristic-PLM exact diagonal entropy wave, density L1 errors on
12, 24, and 48 square grids are `2.09551654e-4`, `7.15524055e-5`, and
`1.61190312e-5`, giving observed orders `1.550234` and `2.150235`.

For characteristic PPM, density errors on 16, 32, and 64 square grids are
`1.44130250e-4`, `3.80399208e-5`, and `9.10858005e-6`, with observed orders
`1.921787` and `2.062216`. At 32 square cells, disabling the transverse
correction increases the characteristic-PPM density error to `4.33952478e-5`.
The exact oblique H2/N2 composition-wave errors are `3.74405262e-5`,
`9.88530404e-6`, and `2.36662586e-6`, corresponding to orders `1.921243` and
`2.062454`.

Both x-normal and y-normal uniform reductions agree with the 1D
characteristic-PLM and characteristic-PPM updates below `3e-12` relative
difference. On the periodic 100-by-4 material-contact regression, bounded
steepening reduces the H2 L1 error from `2.65683522e-4` to `1.66701302e-4`.
The oblique pressure-ratio-three shock remains positive and conservative with
zero pressure overshoot. Enabling flattening changes the conserved solution by
`1.55397791e2` in the pinned mean-absolute signature without increasing total
pressure variation.

The characteristic-PPM 24-square reacting hotspot reaches `2e-6 s` in 59
steps. It retains a maximum conservation error of `2.40e-15`, produces a
`2.82161 Pa` pressure span and `6.91813e-3 m/s` maximum speed, and maintains
roundoff-scale composition closure.


## Qualified dilute-gas molecular transport

The current seven-species transport database stores Lennard--Jones well depth
`epsilon/k_B`, collision diameter `sigma`, geometry, dipole, polarizability, and
rotational-relaxation provenance. The active 0.16.0 subset uses the nonpolar
dilute-gas collision integrals

```text
Omega_mu(T*) = 1.16145/T*^0.14874
             + 0.52487/exp(0.77320 T*)
             + 2.16178/exp(2.43787 T*)

Omega_D(T*)  = 1.06036/T*^0.15610
             + 0.19300/exp(0.47635 T*)
             + 1.03587/exp(1.52996 T*)
             + 1.76474/exp(3.89411 T*)
```

with Chapman--Enskog pure viscosity

```text
mu_k = 2.6693e-6 sqrt(W_k T) / (sigma_k^2 Omega_mu)
```

and binary diffusion

```text
D_kj = 1.8580e-7 T^(3/2) sqrt(1/W_k + 1/W_j)
       / (p_atm sigma_kj^2 Omega_D).
```

Wilke's rule mixes viscosity. Pure conductivity uses a modified Eucken
relation, and the mixture conductivity is the arithmetic/harmonic Mathur mean.
The mass-based mixture-averaged diffusion coefficient is

```text
D_k = (1 - Y_k) / sum_(j != k) X_j / D_kj.
```

Trace regularization prevents singular pure-component denominators without
changing resolved compositions. These formulas are compared against Cantera at
four representative states, but are not labeled full PelePhysics transport
parity because the current implementation omits generated polynomial fits,
polar corrections, and detailed rotational/vibrational conductivity.

## One-dimensional reactive diffusion flux

At each face the adjacent conserved states are converted through the NASA7 EOS.
Transport coefficients are evaluated at arithmetic face temperature, pressure,
density, and composition. For an x-normal face, the viscous stresses are

```text
tau_xx = 4/3 mu du/dx
tau_xy = mu dv/dx
tau_xz = mu dw/dx.
```

The conservative diffusive momentum flux is `-tau`. The energy flux contains
viscous work and Fourier conduction,

```text
F_E,visc+cond = -(tau_xx u + tau_xy v + tau_xz w) - lambda dT/dx.
```

The preliminary ideal-mixture species flux follows the PeleC/PelePhysics form

```text
j_k* = -rho D_k [dX_k/dx + (X_k - Y_k) d(ln p)/dx].
```

A correction velocity is applied through

```text
j_k = j_k* - Y_k sum_j(j_j*),
```

then the final species is assigned the roundoff closure residual so
`sum_k(j_k)=0`. Species enthalpy transport is added to total energy as
`sum_k h_k j_k`. Soret, Dufour, and multicomponent diffusion are excluded.

The diffusion operator is advanced with explicit SSPRK2. Its timestep is
limited by the largest active kinematic viscosity, thermal diffusivity, or
species diffusion coefficient,

```text
dt_transport = C_transport dx^2 / max(4 mu/(3 rho), lambda/(rho cv), D_k),
C_transport <= 0.5.
```

When chemistry and transport are both active, the symmetric composition is

```text
reaction(dt/2) -> transport(dt/2) -> hydro(dt)
               -> transport(dt/2) -> reaction(dt/2).
```

The analytical transverse-shear diffusion wave verifies second-order spatial
and temporal behavior. Separate periodic species and temperature waves verify
smoothing, positivity, species closure, and conservative flux divergence.

## Scope limitations

The code remains serial and uniform-grid. Reactive CFD currently uses only the
seven-species, four-reaction elementary subset, selectable Rusanov/HLLC fluxes,
and qualified frozen-composition PLM/PPM characteristic bases. Molecular
viscosity, Fourier conduction, and mixture-averaged species diffusion are
qualified only in 1D. The code still lacks third-body/falloff chemistry in CFD,
a stiff coupled cell integrator, two-dimensional transport, Soret/multicomponent
transport, a complete mechanism, full general-EOS PeleC Riemann/PPM parity,
complete multidimensional PPM corner tracing, physical boundaries, AMR, MPI,
and accelerators.


## Two-dimensional molecular transport

The 0.17.0 path uses the Newtonian stress tensor, Fourier heat conduction, and
mixture-averaged species fluxes with optional barodiffusion. A correction
velocity enforces zero net diffusive mass flux and species enthalpies contribute
to total-energy transport. SSPRK2 advances the diffusion operator.

## Physical boundary conditions

A no-slip ghost velocity is reflected about the prescribed wall velocity; a
slip ghost reflects only the normal component. Solid-wall inviscid flux is
pressure-only. Adiabatic walls copy temperature and isothermal walls use
\(T_g=2T_w-T_i\). Species diffusive flux is zero at a solid wall. Fixed
inflow uses the configured initial primitive state and outflow uses constant
extrapolation. Periodic boundaries must occur in matched pairs.
