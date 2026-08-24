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

### Hybrid WENO edge reconstruction in AMR

The multilevel AMR path may set `amr_hybrid_weno = .true.` while retaining
`amr_reconstruction = "characteristic_ppm"`. For each primitive five-cell
stencil, the ordinary limited PPM edge pair is replaced by the selected PeleC
formula: WENO5-JS (`amr_weno_scheme = 0`), WENO5-Z (`1`), WENO7-Z (`2`), or
WENO3-Z (`3`). The reconstructed edges still pass through the same
physical-state sanitization, parabolic profile integration, and characteristic
projection described above. As in PeleC's hybrid branch, shock flattening is
part of the ordinary PPM alternative rather than an extra blend on WENO edges.

This option uses the existing four coarse/fine ghost layers, midpoint
parent-time interpolation, recursive SSPRK3 stages, and effective reflux flux.
It is qualified only for the one-dimensional multilevel reactive path.
Regular-grid 2D use, multiple AMR patches, and distributed patch ownership are
not included in this milestone.

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

## Pressure-dependent reaction rates

Third-body rates use the efficiency-weighted collider concentration. Falloff rates combine low- and high-pressure limits with the reduced-pressure factor and optional Troe broadening. The full CFD path uses adaptive implicit backward-Euler trials with step doubling and Richardson extrapolation.

## Distributed one-dimensional operators

The MPI path uses uneven contiguous blocks and one ghost cell on either side of
each local state. Periodic nonblocking halo exchange supplies face states for
hydro and molecular transport. `MPI_Allreduce` synchronizes the hyperbolic and
parabolic timestep limits, operator success flags, and conservation diagnostics.

The coupled update is transactional: chemistry, transport, hydro, transport,
and chemistry act on a trial copy. If any local implicit solve or state recovery
fails, the failure is reduced globally, the trial is discarded on every rank,
and the complete Strang interval is retried at half size. This keeps rank-count
changes from producing divergent accept/reject histories.

## Scope limitations

The current distributed path is one-dimensional and uniform-grid. Serial 2D
reactive flow includes molecular transport and physical boundaries, while the
AMR layer provides arbitrary-depth 1D reactive PCM/PLM/PPM/WENO advancement,
molecular transport, synchronization, tagging, overlap-preserving single-patch
dynamic regridding, outflow-boundary refinement, and composite output. MPI 2D,
dynamic multipatch application integration and depth, one-sided periodic-seam
refinement, embedded boundaries, LES, particles/spray, and accelerators remain
outside the implemented scope. Fixed two-level separated multipatch WENO hydro,
chemistry, molecular transport, and conservative regrid transfer are qualified.
The
transport model still
excludes Soret, Dufour,
multicomponent Stefan--Maxwell diffusion, polar corrections, and bulk viscosity.
The characteristic hydro basis remains a qualified frozen-composition
approximation rather than complete PeleC/PelePhysics general-EOS Riemann and
multidimensional PPM corner-tracing parity.

## Static two-level AMR transfer and synchronization

For a coarse cell average `U_i` and refinement ratio `r`, the fine child center
has normalized coarse-cell offset

```text
xi_j = (j - 1/2) / r - 1/2,  j = 1 ... r.
```

The first AMR prolongation uses an MC-limited coarse slope `s_i`:

```text
U_f(i,j) = U_i + s_i xi_j.
```

Because the child offsets sum to zero, their arithmetic mean is exactly `U_i`.
Restriction therefore uses the matching volume average. The fine patch is
strictly interior in this slice so every prolonged coarse cell has both slope
neighbors and the patch has two explicit coarse/fine interfaces.

For each interface the flux register stores the time-integrated mismatch

```text
delta I = sum_fine(dt_f F_f) - dt_c F_c.
```

At each existing coarse/fine interface, reflux subtracts `delta I / dx_c` from
the uncovered coarse cell on the left or adds it to the uncovered coarse cell
on the right. A side coincident with a physical boundary has no uncovered
parent neighbor and receives no reflux. Covered coarse cells are separately
replaced by restricted fine averages. A composite integral counts uncovered
coarse volumes and fine volumes exactly once, providing the conservation gate.

## Arbitrary-depth nested AMR hierarchy

An arbitrary-depth hierarchy is an allocatable chain of adjacent-level
relations. For interface `ell -> ell+1`, the relation records its own integer
ratio `r_ell`, parent-local patch bounds, physical parent bounds, and parent and
child spacing. The child becomes the complete parent field of the next
relation, so no fixed maximum level count or common refinement ratio is built
into the data model.

The number of substeps and step size relative to a root interval are

```text
N_0 = 1,
N_ell = product(m=0...ell-1, r_m),
dt_ell = dt_0 / N_ell.
```

Conservative prolongation proceeds from root to deepest level. Restriction,
reflux, and average-down proceed in the opposite direction. For each relation,
reflux first corrects its one or two existing coarse/fine sides, then
average-down replaces the covered parent cells. Processing the deepest
interface first ensures the state transferred to its parent already contains
all finer-level information.

The multilevel composite integral sums the uncovered portion of every level
that owns a child and the complete deepest field. A four-level gate with ratios
`2`, `3`, and `2` verifies cumulative geometry, 12-way root-relative
subcycling, restriction identities, and simultaneous conservation correction
at all three interfaces.

## Arbitrary-depth reactive AMR advance

`amr_multilevel_reactive_1d_mod` owns one conserved field, temperature field,
and ghost pair for every runtime level. Given cumulative refinement
`R_ell = product(r_0...r_(ell-1))`, the root timestep satisfies

```text
dt_0 <= min_ell(R_ell * dt_hydro,ell),
dt_0 <= min_ell(R_ell^2 * dt_transport,ell).
```

Hydro advances a parent over its interval, saves its beginning and provisional
end states, and advances the child `r_ell` times with interpolated parent ghost
states. Each child call recursively completes all deeper work before returning
its time-integrated outer fluxes. The caller then refluxes and averages down its
own parent/child relation. Molecular transport uses the same recursion with
`r_ell^2` child substeps and the actual parent/child center distance at outer
child faces.

Cell-local chemistry advances every level over the same physical half interval,
followed by deepest-to-root average-down. The complete accepted update is

```text
R(dt/2) -> T_recursive(dt/2) -> H_recursive(dt)
        -> T_recursive(dt/2) -> R(dt/2).
```

The solution is deep-copied before this composition, so any failed chemistry,
transport, hydro, EOS recovery, reflux, or ghost fill restores every level,
the root time, and the step counter. A three-level hotspot gate exercises PLM,
molecular transport, chemistry, recursive synchronization, composite
conservation, positivity, and species closure.

### Tag-driven multilevel regridding and output

Starting from the root, the runtime application evaluates the existing
normalized-gradient criterion on each parent. If tags exist, it builds one
patch, conservatively prolongs that child, and repeats until no tags remain or
`amr_max_levels` is reached. Internal coarse/fine sides retain the parent-cell
buffer required by the four wide ghost layers and limited prolongation. An
outflow side coincident with the global physical boundary may instead extend
to that edge at every nested level.

At a physical outflow side, the fine face-adjacent and wide PPM/WENO ghosts use
constant extrapolation from the fine boundary cell. The flux register ignores
that side because no uncovered parent neighbor exists; reflux is applied only
at the opposite coarse/fine interface. A periodic fine patch may use fine-level
wrap only when it covers the complete parent domain. One-sided periodic-seam
refinement is rejected.

Before changing a hierarchy, every child is averaged down deepest-to-root. The
new nested chain is then planned from this synchronized root and every child is
conservatively prolonged. Thus the complete composite integral is invariant
under creation, movement, resizing, depth reduction, and depth growth. If the
new bounds and ratios are identical, the existing hierarchy is retained
without reconstruction. If they differ, every common old/new fine level is
intersected in physical coordinates. When its spacing agrees and both overlap
edges lie on cell boundaries, conserved state and temperature are copied
exactly into the rebuilt patch. Processing all levels and then averaging down
deepest-to-root makes retained deepest data authoritative in covered parents.
Levels whose spacing changed retain their conservative prolongation instead.

Composite output follows the hierarchy recursively:

```text
write uncovered parent left -> write child composite -> write parent right.
```

This traversal emits monotonically ordered cell centers and counts each domain
volume exactly once at its finest active representation.

## One-dimensional tagging and dynamic regridding

For a selected component `q`, cell `i` uses the largest adjacent jump,

```text
jump_i = max(|q_i - q_(i-1)|, |q_(i+1) - q_i|),
scale_i = max(scale_floor, |q_(i-1)|, |q_i|, |q_(i+1)|).
```

A cell is tagged when the jump is above the absolute threshold and
`jump_i / scale_i` is at least the relative threshold. The one-patch planner
takes the bounding interval of all tags, adds a configured coarse-cell buffer,
and expands deterministically to the requested minimum width.

Before replacing an existing patch, its fine values are averaged down. The new
patch is then conservatively prolonged from that synchronized coarse state. If
the refinement ratio is unchanged, fine cells in the geometric overlap are
copied from the old patch after prolongation. Thus leaving regions preserve
their fine averages, entering regions preserve their coarse averages, and the
overlap preserves its complete fine representation. This sequence makes the
composite integral invariant under patch creation, movement, resizing, and
removal.

## Reactive two-level AMR advance

The reactive AMR application uses PCM or limited primitive PLM face states and
the qualified general-EOS Rusanov or HLLC flux. If the fine-level explicit
limit is
`dt_f,max`, the coarse interval is limited by

```text
dt = min(dt_coarse,max, r * dt_f,max),
dt_f = dt / r.
```

The coarse level advances once and the fine patch advances `r` times. Fine
ghost states at substep `m` interpolate the adjacent uncovered coarse state
between the beginning and end of the coarse hydro interval with
`alpha = m/r`. Coarse and fine fluxes at both patch interfaces are accumulated
over their actual time steps, then the existing reflux correction is applied.

Chemistry is composed symmetrically around the complete hydro interval,

```text
R(dt/2) -> A_coarse(dt) + [A_fine(dt/r)]^r -> reflux/average-down
        -> R(dt/2) -> average-down.
```

The cell-local constant-volume reactor conserves density, momentum, and total
energy. Reflux accounts for coarse/fine advective mismatch, and average-down
counts the fine solution in covered volumes, so the composite mass, momentum,
and energy conservation argument remains valid for the coupled reactive step.

### Limited PLM and flux-consistent SSPRK2

For AMR PLM, each primitive component uses the configured limited slope. The
left and right cell-edge values are

```text
q_i,L = q_i - s_i/2,
q_i,R = q_i + s_i/2.
```

Density and pressure fall back to the cell center if reconstruction violates a
floor. Species face values are clipped at zero and normalized to unit sum before
the general-EOS primitive-to-conserved conversion.

With `L(U)` denoting the flux divergence, the level update is SSPRK2,

```text
U(1)   = U(n) + dt L(U(n)),
U(n+1) = 1/2 U(n) + 1/2 [U(1) + dt L(U(1))].
```

This update is exactly conservative with the effective face flux
`F_eff = [F(U(n)) + F(U(1))]/2`. The AMR flux register therefore accumulates
`F_eff`, rather than only one stage, at both coarse/fine interfaces. PLM fine
substeps use the time-interpolated coarse ghost state at each substep midpoint.

### Molecular transport and diffusive synchronization

The AMR transport operator uses the uniform-grid 1D diffusive face kernel for
Newtonian stress, Fourier conduction, mixture-averaged species diffusion,
optional barodiffusion, correction velocity, and species-enthalpy transport.
It is composed symmetrically with reaction and advection:

```text
R(dt/2) -> T(dt/2) -> A(dt) -> T(dt/2) -> R(dt/2).
```

For a parabolic fine-level stability limit `dt_T,f,max`, the coarse interval is
limited by

```text
dt <= min(dt_T,c,max, r^2 dt_T,f,max).
```

Each transport half interval advances the coarse level once with SSPRK2 and
the fine level in `r^2` SSPRK2 substeps. Fine ghosts use coarse states
interpolated to each substep midpoint. Interior fine faces use `dx_f`; a patch
interface uses the actual coarse-center to fine-center distance
`(dx_c + dx_f)/2` when forming its gradient.

The effective diffusive flux is the arithmetic mean of the two SSPRK2 stage
fluxes. Its time integral is accumulated at each coarse/fine interface, and a
dedicated reflux correction is followed by covered-cell average-down. Because
the species correction velocity makes the total diffusive mass flux zero, the
same synchronization conserves total mass and each periodic species mass while
the energy flux includes conduction, viscous work, and transported species
enthalpy.

## Separated multipatch AMR

A patch set stores ordered coarse index intervals
`[lo_p, hi_p]` with at least one uncovered parent cell between consecutive
patches. This separation makes every nonphysical patch side a true coarse/fine
interface. Adjacent tag candidates are coalesced before hierarchy creation;
same-level fine/fine ghost exchange is therefore not yet required.

The composite integral is

```text
I = dx_c * sum(uncovered parent cells U_c)
  + dx_f * sum_p sum(all cells of patch p U_f,p).
```

Each covered parent interval is subtracted exactly once. Prolongation,
restriction, flux-register correction, and average-down are applied per patch
inside a transactional set-wide operation. At regrid, all old patches first
average down to the parent. Every new patch is conservatively prolonged, then
each old/new pair copies its same-ratio fine intersection. This permits patch
movement, split/repartition, creation, and removal without losing aligned fine
structure or changing the composite integral.

For the qualified fixed two-level reactive path, the parent hydro advances
once over `dt`. Every patch advances independently for `r` steps of `dt/r`,
using the same parent start/end states for time-interpolated narrow and
four-layer PPM/WENO ghosts. One flux register is accumulated per patch. All
registers are refluxed and every covered region is averaged down before the
parent temperature is recovered. Each molecular-transport half interval
advances the parent once and every patch for `r^2` steps, accumulating the
mean SSPRK2 diffusive face flux in its own register. Chemistry advances all
parent and patch cells for the same physical half interval and then averages
down the patch set. The full transactional composition is
`R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)`. The root timestep is limited by `r`
times every fine hyperbolic limit and `r^2` times every fine parabolic limit.
At each configured regrid interval, the public two-level multipatch driver tags
the synchronized parent, clusters tags separated by no more than
`amr_maximum_patch_gap_cells`, and expands each cluster by the configured buffer
and minimum width. The old set first averages down; the new set is
conservatively prolonged; every equal-resolution old/new fine intersection is
then copied exactly. An empty tag collection removes all patches. Ordered
composite output alternates uncovered parent intervals and fine patch cells so
the domain is covered exactly once.

For a static arbitrary-depth patch tree, relation `l` stores a separated child
set for every patch at level `l-1`. Child offsets flatten those parent-local
sets into one deterministic level order. Conservative prolongation walks from
the root toward the leaves. Average-down walks from the deepest relation to the
root, so a changed leaf is restricted through every ancestor. The composite
integral begins with the root, then for every child subtracts its covered
parent interval and adds its fine-cell integral. This replacement formula is
valid for branching trees and mixed per-level refinement ratios. Reactive
flux registers mirror the geometry: each relation contains one register array
per parent and one register per local child. Synchronization applies every
deepest relation's reflux before restricting that relation into its parents,
then continues toward the root. Fields and registers roll back together if any
parent-set operation fails.

For reactive patch-tree hydro, a recursive call advances one parent patch for
`dt_l` and records the time-integrated coarse flux at every local child
boundary. Each child then advances recursively for `r_l` intervals of
`dt_l/r_l`. Child ghost cells interpolate between the parent state before and
after its advance. The returned fine boundary-flux integrals accumulate in the
matching child register. Once all child substeps finish, the parent-owned set
is refluxed and averaged down. Because descendants are synchronized before a
child returns, this ordering propagates every leaf correction through all
ancestors. The stable root interval is the minimum of each patch's local CFL
limit multiplied by its cumulative refinement ratio.

The complete reactive-tree interval uses
`R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)`. Each reaction operator advances every patch for the
same physical half interval because chemistry is cell-local rather than a
mesh-stability-limited flux operator. It then averages the changed fields from
the deepest leaves to the root, recovers every parent temperature, and refreshes
all ghosts. A transaction around the full composition prevents a failed second
reaction half-step from retaining the already accepted hydro state.

Each transport recursion advances its parent for `dt_l`, then advances each
child for `r_l^2` intervals of `dt_l/r_l^2`. The parent start and end states
provide midpoint-in-time coarse/fine ghosts. Coarse and fine diffusive boundary
fluxes accumulate in the same nested register shape as hydro, followed by
per-parent reflux and average-down. A fine transport limit is converted to a
root interval by the square of its cumulative refinement ratio.

The qualified reactive-tree gate uses PCM, elementary chemistry, molecular
transport, and strictly interior separated patches.

For a runtime patch-tree rebuild, the old tree is first averaged down so its
root represents the complete composite solution. The new branching plan is
conservatively prolonged from that root. At each common level with equal `dx`,
the physical intersection of every old/new patch pair is converted to aligned
cell offsets and copied exactly. This coordinate rule survives patch movement,
repartition, and changed parent ownership. The rebuilt tree is then averaged
down deepest-to-root.

For automatic tree planning, the synchronized root is the canonical input.
Each prospective parent is tagged with the configured normalized-gradient and
absolute thresholds. Its disconnected tag groups are buffered, expanded to a
minimum width, and assigned that parent's flattened index. Children are
conservatively prolonged before the same operation repeats at the next
relation. Recursion stops independently on untagged branches and globally when
no children remain or `amr_max_levels` is reached. Physical-boundary children
remain a future integration.

Adjacent children of one parent remain separate owners. Before each fine
substep, a child first receives its time-interpolated parent ghosts. Any ghost
fine index covered by a sibling is then replaced exactly by that sibling's
interior state and temperature; the same lookup fills up to four PPM/WENO
layers across a chain of adjacent boxes. Each patch advances locally and
returns its time-integrated boundary flux. The two values at a shared
fine/fine face are replaced by their arithmetic mean, and conservative
boundary-cell corrections make both patches use that single owned flux. The
corresponding register sides are zeroed before reflux because neither is a
coarse/fine interface. Hydro and molecular transport use the same interface
ownership rule.

## MPI AMR patch distribution bridge

For each patch in level/flattened-patch order, `0.49.0` assigns the patch to
the rank with the smallest accumulated owned cell count. The root patch enters
the same schedule, and a tie selects the lowest rank, making the owner map
deterministic once the hierarchy and communicator size are fixed. Before
assignment, communicator-wide minimum/maximum reductions require identical
base cells, relation ratios, parent/child topology, child bounds, and physical
root extent on every rank.

Patch values use an owner-authoritative replicated bridge:

```text
owner patch values --MPI_Bcast--> every rank's patch replica
```

For an adjacent sibling pair, the owner of the left patch broadcasts its last
`g` cells and the owner of the right patch broadcasts its first `g` cells.
Those buffers become the opposite sibling's ordered halo layers, with layer
one nearest the shared face. All ranks enter broadcasts in the same hierarchy
order, so the schedule is independent of the local owner set.

The `0.50.0` chemistry operator advances each patch on its owner only. After
each local reactor call, `MPI_Allreduce(MPI_LAND)` accepts or rejects the
result collectively. An accepted owner broadcasts state, temperature, and
narrow/wide ghost storage; a rejection restores the owner-synchronized backup
on every rank. After the last patch, the replicated tree averages down
deepest-to-root, recovers thermodynamic temperatures, and rebuilds ghosts.
Thus each global chemistry interval contains exactly one reactor call per
patch independent of communicator size.

The present bridge still stores all patch interiors on every rank and does not
claim sparse storage, owner-only hydro/transport, distributed shared-flux
correction, or regrid migration.
