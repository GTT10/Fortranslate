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
\(T_g=2T_w-T_i\). Species diffusive flux is zero at a solid wall by default.
A prescribed wall may instead provide wall-to-gas mass fluxes \(J_k^w\) with

```text
sum_k J_k^w = 0,
F_rhoY,k(lower) =  J_k^w,
F_rhoY,k(upper) = -J_k^w,
F_E,species = sum_k h_k F_rhoY,k.
```

The zero-sum constraint preserves the impermeable total-mass and momentum
wall contract while permitting species conversion. The species positivity
limiter scales the complete vector and its enthalpy flux together. Fixed
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

The `0.51.0` hydro operator walks the recursive patch-tree schedule
collectively. On each patch interval, only its owner calls the same one-patch
PCM/PLM/PPM finite-volume kernel used by serial recursion. After logical-AND
acceptance, the owner broadcasts the interval-start state, all face fluxes,
the accepted reactive patch, and the updated level counter. All replicas can
therefore construct the same interval-end state and boundary flux integrals.
They fill time-interpolated parent ghosts, overwrite adjacent sibling halos,
advance child intervals in deterministic order, reconcile each shared
fine/fine time-integrated flux, accumulate coarse/fine registers, reflux, and
average down. A rejected patch or synchronization step restores the complete
owner-synchronized solution and reports zero accepted local calls.

At `0.51.0`, the bridge still stores all patch interiors on every rank. Shared
flux reconciliation and hierarchy synchronization are deterministic replicated
operations after owner broadcasts, rather than a sparse stage-synchronous
exchange. That milestone does not claim owner-only molecular transport,
sparse storage, regrid migration, or scalable point-to-point communication.

The `0.52.0` transport operator applies the same ownership rule to the shared
one-patch SSPRK2 transport kernel. Viscosity, Fourier conduction,
mixture-averaged species diffusion, barodiffusion, correction velocity, and
species-enthalpy flux remain coupled exactly as in serial AMR. Each refinement
relation performs `r^2` child transport intervals, and owner-weighted call
counts therefore use the square of every cumulative refinement ratio. Owners
broadcast the interval-start state and the SSPRK2 effective face flux before
replicas apply adjacent shared-flux correction, diffusive reflux, and
average-down. Failure at any depth restores the synchronized pre-call tree.

The transport qualification remains owner-authoritative and replicated. It is
not sparse stage communication and is exposed as a separate distributed
operator; a single outer chemistry--transport--hydro transaction, regrid
migration, and scalable point-to-point schedules remain pending.

The `0.53.0` distributed reactive step is

```text
owner chemistry(dt/2)
owner molecular transport(dt/2)
owner recursive hydro(dt)
owner molecular transport(dt/2)
owner chemistry(dt/2)
```

The complete synchronized solution is copied before the first operator. Each
inner operator reaches communicator-wide acceptance independently, but only
the outer call commits the sequence. A missing transport database is rejected
before mutation. If, for example, chemistry and the first transport interval
succeed but hydro rejects its reconstruction, every rank restores state,
temperature, ghosts, time, step, advance counters, and regrid statistics to
the outer backup. Returned owner-call counts are zero for a rejected sequence.

Successful call accounting contains two owner visits per patch for chemistry,
one cumulative `r` schedule for hydro, and two cumulative `r^2` schedules for
transport. The current implementation composes owner-authoritative replicated
operators; it does not yet provide sparse storage or patch migration.

## Embedded-boundary geometric fractions

Let nodal `phi > 0` denote fluid and `phi <= 0` denote solid/interface. Each
Cartesian cell is divided into the triangles `(00,10,11)` and `(00,11,01)`.
Within a triangle, `phi` is affine. Edges that change sign are intersected by
linear interpolation, the positive polygon is clipped, and its shoelace area
contributes to the cell volume fraction. The two triangle areas sum to a value
in `[0,1]`.

The open fraction of a Cartesian face is the positive length of its linearly
interpolated endpoint values. If the positive normalized interval is
`[s0,s1]`, the open fraction is `a=s1-s0` and the stored tangential centroid
offset is `c=(s0+s1)/2-1/2`. Thus `-0.5 <= c <= 0.5`; x-face offsets are
normalized by `dy` and y-face offsets by `dx`, matching the AMReX face-centroid
convention. Full and closed faces store zero offset. Fractions within the
roundoff classification tolerance become covered or regular; all intermediate
values are cut cells.
This construction reproduces any planar interface exactly on the mesh and
gives second-order integrated-area convergence for the circular level-set
gate. Each triangle's zero contour supplies a physical segment, centroid, and
normalized `grad(phi)`. A cut cell combines its two segments by physical-length
weighting; coincident diagonal segments are counted once. The resulting unit
normal points from solid toward fluid. A separate integrated normal preserves
the segment-vector sum used by the pressure force. Cartesian fluxes are
multiplied by their shared open fractions, combined with that wall force, and
divided by the fluid volume.

The first-order small-cell path follows FluxRedist. For cut-cell volume fraction
`kappa`, it forms `Rnc` by volume-weighting the conservative right-hand side
over the cell and its positive-aperture face neighbors. The local stable update
is `kappa*Rc + (1-kappa)*Rnc`. The removed extensive update
`kappa*(1-kappa)*(Rc-Rnc)` is divided by the sum of neighbor volume fractions
and added to each neighbor. Thus the volume-weighted domain update is unchanged
while the cut-cell contribution no longer contains an unscaled `1/kappa`
factor. The reactive forward update is committed only after every active cell
passes conserved-to-primitive EOS recovery.

The alternative weighted StateRedist path operates on the provisional state,
not on its right-hand side. A cell with `0 < kappa < 0.5` selects the dominant
direction of the aperture-difference normal. If the accumulated fluid volume
is below `0.5`, or the two normal components have equal magnitude, the
orthogonal neighbor and then the connecting diagonal are included. This is the
two-dimensional AMReX `MakeITracker` construction with at most three
neighbors and no periodic wrapping.

Let `M_i` be the selected neighbors of cell `i`, and let `N_j` count how many
neighborhoods contain cell `j`, including its own. For a merging cell,
`beta_i = (0.5-kappa_i)/sum_(j in M_i) kappa_j`. Its self partition is
`alpha_i = 1-sum_(d: i in M_d) beta_d/N_i`. The weighted neighborhood volume
and state are
`Vhat_i = alpha_i*kappa_i + sum_(j in M_i) beta_i*kappa_j/N_j` and
`Qhat_i = (alpha_i*kappa_i*U_i + sum_(j in M_i)
beta_i*kappa_j*U_j/N_j)/Vhat_i`. The default `max_order=0` redistribution
returns `U'_j = alpha_j*Qhat_j + sum_(i: j in M_i) beta_i*Qhat_i/N_j`.
The partition makes the volume-weighted sum of `U'` equal that of `U`, while a
uniform state gives every `Qhat` and every output cell the same value.

For `max_order=2`, the same partition locates `Qhat_i` at
`Xhat_i = (alpha_i*kappa_i*c_i + sum_(j in M_i)
beta_i*kappa_j*(j-i+c_j)/N_j)/Vhat_i`, where `c` is the normalized
fluid-volume centroid offset. A least-squares linear fit of neighboring
`Qhat` values uses the connected 3-by-3 stencil and grows to active cells in a
5-by-5 stencil when the normal matrix is rank deficient. Pairwise predictions
at neighboring `Xhat` locations receive the AMReX centroid limiter. A final
common slope scale bounds every self and merge-recipient evaluation by the
active input component range. Because `Xhat_i` is the weighted centroid of the
same recipient partition, the linear corrections have zero volume-weighted
first moment; conservation and affine-state reproduction therefore hold
together.

The reactive StateRedist advance builds `U = U_old + dt*R`, redistributes all
conserved components together, and commits state and recovered temperature
only if every active cell passes the general-EOS conversion. Fourth-order
StateRedist slopes, periodic/ghost-cell neighborhoods, and multilevel
redistribution are not yet claimed.

## Static two-level EB average-down

One fine rectangle may be aligned with a coarse EB level at integer refinement
ratio `r`. Before transfer, the hierarchy verifies matching physical bounds and
spacing and requires each covered parent volume fraction to equal the arithmetic
average of its `r^2` child fractions. For a parent cell with positive total fine
fluid measure, component `n` is restricted as

`U_c(n) = sum_f(kappa_f U_f(n)) / sum_f(kappa_f)`.

If every child is covered, the first child value supplies the generic fallback,
matching the AMReX EB average-down kernel. The reactive wrapper instead retains
the original covered-parent state and temperature, while every active parent
must pass conserved-to-primitive EOS recovery before the complete arrays are
committed.

The composite extensive diagnostic excludes every coarse cell geometrically
covered by the fine rectangle and adds `kappa_f U_f dx_f dy_f` over the fine
level. Restricting the fine patch and then integrating the resulting coarse
level therefore reproduces the same composite integral to roundoff. This
foundation does not advance either level and does not yet provide EB
prolongation, ghost fill, subcycling, flux-register reflux, dynamic regridding,
multiple patches, or distributed ownership.

## EB coarse/fine flux register and re-reflux

The flux register accumulates coarse and fine time intervals independently.
For a low-x interface, the exterior coarse cell receives the raw state
correction

`D = sum_c(dt_c a_c F_c)/(kappa_c dx_c)
     - sum_f(dt_f a_f F_f dy_f)/(kappa_c dx_c dy_c)`.

High-x reverses both signs; the y interfaces use the corresponding `dy_c` and
fine `dx_f` measures. This form supports fine subcycling without assuming that
coarse and fine call counts match. When flux and aperture data match physically,
one coarse interval and all fine subintervals cancel to roundoff.

A regular exterior cell receives `D` directly. For a cut cell, define the raw
extensive correction `dm = kappa_c D`. Re-reflux adds `dm` to that cell's state,
then adds `dm(1-kappa_c)/sum_n(kappa_n)` to every connected neighbor in its
3-by-3 stencil. The resulting fluid-volume-weighted correction is exactly
`kappa_c D`, while the small cell itself receives only its stable `kappa_c`
share. A neighbor geometrically covered by the fine rectangle is represented
by all of its fine children; adding the same state increment to each child
preserves the parent measure because their volume fractions restrict to the
coarse value.

The generic operation commits both state arrays and resets the register only
after all corrections are finite. The reactive wrapper operates on a copied
register, restores covered-cell data, recovers active coarse and fine
temperatures, and commits the two levels only after every EOS conversion
succeeds. This is interface synchronization only; it does not yet compose EB
level advancement, ghost fill, prolongation, or regridding.

## Static two-level reactive EB hydrodynamic advance

Piecewise-constant prolongation initializes child `f` from its parent `c` as
`U_f=U_c`. Covered child data is retained as injected; every active child is
accepted only after conserved-to-primitive EOS recovery supplies a valid
temperature.

One coarse interval uses the sequence

`A_c(dt) -> [G_c(alpha_m), A_f(dt/r)] for m=1..r -> reflux -> average-down`,

where `A` is the existing centroid-flux and weighted-StateRedist EB hydro
operator. For PCM, `alpha_m=(m-1)/r`; for the time-centered characteristic-PLM
trace, `alpha_m=(m-1/2)/r`. At every open fine-patch boundary face, `G_c`
forms the adjacent exterior state

`U_ext(alpha)=(1-alpha) U_c^n + alpha U_c^(n+1)`

and recovers its temperature with the mixture EOS. Spatial transfer is
piecewise constant: the `r` fine faces belonging to one coarse face sample the
same adjacent coarse cell. Closed faces never consume their placeholder
exterior values.

The coarse step contributes its actual centroid fluxes for `dt`; each fine
substep contributes its own actual centroid fluxes for `dt/r`. The established
EB re-reflux then corrects the coarse/fine mismatch, and reactive EB
average-down replaces the refined parent region. Caller outputs remain equal
to their inputs unless every level advance, accumulation, EOS recovery,
re-reflux, and restriction succeeds. This qualified composition covers one
strictly internal aligned rectangle and hydrodynamics only; it does not yet
include coarse spatial slopes, physical-domain-touching fine patches,
chemistry, molecular transport, regridding, multiple patches, deeper levels,
or MPI ownership.

## Runnable static EB AMR timestep selection

For each accepted hierarchy step, the driver evaluates the existing active-cell
hyperbolic CFL limit independently on both levels. If `dt_c` and `dt_f` are
those limits and the refinement ratio is `r`, the public coarse interval is

`dt = min(dt_c, r dt_f, t_final-t)`.

The fine operator therefore receives `dt/r <= dt_f` on every one of its `r`
substeps. Stability is recomputed from the synchronized hierarchy after each
accepted interval. A maximum-step bound rejects an unfinished run, and no
output hierarchy is reported as successful unless the final time is reached
with a positive finite minimum accepted timestep.

The input defines one strictly internal rectangle with inclusive one-based
coarse indices. Its physical bounds follow exactly from the root spacing. The
same analytic plane or circle level set is sampled on the root nodes and the
refined patch nodes before the existing geometry-measure compatibility gate is
applied. Initial PCM prolongation, every time interval, and final composite
diagnostics reuse the qualified two-level transactions. Separate CSV files
retain the complete synchronized parent and child layouts for inspection.

The complete EB hydro update accepts either `pcm` or `characteristic_plm`.
PCM supplies the cell primitive state unchanged. Characteristic PLM reuses the
regular reactive frozen-composition acoustic projection and MUSCL-Hancock
normal trace. A cell receives a limited normal slope only when both stencil
neighbors are active; otherwise that directional slope is exactly zero. The
density, pressure, and species bounds are applied before the trace, and each
traced composition is clipped and normalized before NASA7 conversion.

The selected reactive Riemann solver first produces `F_center` on every
positive-aperture Cartesian face. A closed face is exactly zero. At a
nonperiodic domain face, the exterior state is the adjacent cell-center state,
giving the qualified zero-gradient boundary. In two dimensions, a partial
face with normalized tangential centroid offset `c` uses

`F_centroid = (1-|c|) F_center(j) + |c| F_center(j+sign(c))`.

The x/y roles are exchanged for a y-face. Interpolation is used only when the
selected tangential face is in-domain and open; otherwise the local
face-center flux is retained so no covered-face value enters the stencil. The
open fraction is applied separately by the conservative divergence, avoiding
double area weighting. The integrated slip-wall pressure force completes `R`,
then weighted StateRedist transactionally advances `U_old + dt*R`.

This is a qualified normal characteristic-PLM plus face-centroid interpolation
path. It does not include PeleC's unsplit EB transverse predictor, PPM, or
fourth-order StateRedist slopes.

The runnable EB driver uses the regular active-cell hyperbolic rate
`(|u|+c)/dx + (|v|+c)/dy` and takes `dt = CFL/max(rate)`, ignoring covered
cells. No inverse-volume-fraction factor is added to this bound because the
weighted StateRedist update removes the explicit small-cell stiffness. The
last step is clipped to the configured final time. Conserved diagnostics use
`sum(kappa*U)*dx*dy`; extrema and EOS checks visit active cells only.

With chemistry enabled, one EB interval is
`reaction(dt/2) -> EB hydro(dt) -> reaction(dt/2)`. The cell-local constant
volume reactor receives a logical active mask, so covered cells are never
converted or integrated. Chemistry operates on candidate state and temperature
arrays and commits only after every active cell succeeds. The outer EB Strang
operator retains the original inputs until both reaction halves and the hydro
transaction complete, so a later Riemann or EOS failure also rolls back the
first reaction half.

## Temperature-tagged EB AMR patch replacement

For each active root cell at least one cell away from the physical boundary,
the regrid indicator is the largest temperature jump to an active Cartesian
neighbor,

`g_ij = max_neighbor |T_neighbor - T_ij|`.

A cell is tagged when `g_ij` exceeds the configured absolute floor and
`g_ij / max(T_ij, T_neighbor, T_floor)` meets the relative threshold. Covered
cells never tag and never enter a neighbor jump. The planner takes one bounding
box around all tags, grows it by the requested buffer, enforces minimum x/y
extents, and clamps it to the strictly internal root region required by the
coarse exterior-state fill.

Replacing old patch `P_old` by `P_new` uses

`average-down(P_old) -> PCM(P_new) -> retain(P_old intersect P_new)`.

The first operation transfers regions that lose refinement into the root. PCM
therefore preserves the root integral in newly refined cells. Exact copying by
global fine index restores all same-resolution overlap after PCM, avoiding
unnecessary diffusion. State and temperature arrays on both levels commit only
after the new patch geometry is valid, every candidate is finite, and EOS
recovery succeeds on every active new fine cell.

Fine-patch removal is optional. On an empty plan, the lifecycle transaction
average-downs the complete child into a candidate root, recovers all affected
active temperatures, and only then commits the root and releases the fine state,
temperature, geometry, and patch metadata. While no fine patch exists, the
coarse timestep is the root active-cell CFL limit and the existing single-level
reactive EB hydro operator advances the state. If tags later return, the driver
builds a new aligned geometry and uses PCM prolongation from the synchronized
root before publishing the child. Conserved diagnostics likewise select either
the two-level composite integral or the root volume-weighted integral. The
default policy continues to retain an untagged patch for backward compatibility.

## Reactive chemistry on the EB AMR lifecycle

For an active two-level hierarchy, one accepted coarse interval is

`R_c,f(dt/2) -> H_EB-AMR(dt) -> R_c,f(dt/2) -> average-down`,

where `R_c,f` applies the cell-local constant-volume reactor independently to
active coarse and fine cells, and `H_EB-AMR` is the existing coarse step, `r`
fine substeps, EB re-reflux, and reactive average-down transaction. Each level
constructs its chemistry mask from the EB cell types, so covered state and
temperature are never passed to the reactor. The final average-down restores
the coarse representation beneath the reacted fine patch.

All four output arrays initially contain the caller's input. Reaction and hydro
operate on private candidates, and those outputs are replaced only after both
reaction halves, EB hydro synchronization, EOS recovery, and final restriction
succeed. Thus a failure after the first reaction half, including an invalid
Riemann solver, rolls back the complete coarse/fine state and temperature. If
the lifecycle has no fine patch, the driver instead applies the already
qualified single-level `reaction-hydro-reaction` EB operator to the root.

This composition contains no molecular transport. Inputs requesting it are
rejected before timestep selection or state mutation.

## Reactive EB AMR checkpoint transaction

The serial checkpoint is a versioned formatted stream containing a magic
header, schema number, species names, configuration signature, actual lifecycle
topology, run metadata, root payload, optional fine payload, and terminal
marker. The configuration signature fixes the mesh and embedded-boundary
definition, refinement ratio, chemistry model and tolerances, Riemann and
reconstruction choices, redistribution settings, and dynamic-regrid policy.
Final time, maximum step count, output paths, and checkpoint scheduling are not
part of the signature and may change for continuation.

Restart first validates the header, schema, species ordering, signature, and
metadata. It reconstructs coarse and fine EB metrics from the input and stored
actual patch bounds, then reads state into private arrays. Conserved state is
authoritative: every active temperature is recovered through the general EOS,
while finite covered data are retained. Only a complete terminal marker and
successful EOS recovery publish the candidate hierarchy. Any parse,
compatibility, geometry, or EOS failure returns invalid, unallocated outputs.

An accepted step writes only after any scheduled regrid transaction. Thus a
checkpoint can encode either an active fine rectangle or a root-only lifecycle
state. Optional stop-after-write exits without forcing time to the requested
final value. On restart, stored step and regrid counters preserve cadence and
the stored minimum timestep preserves the cumulative diagnostic.

## Two-level reactive EB patch sets

The patch-set planner scans tagged root cells in deterministic index order and
clusters cells connected within a Chebyshev reach of
`maximum_patch_gap_cells + 1`. Each component receives its own buffered,
minimum-size bounding rectangle. Rectangles that would leave fewer than two
coarse cells of separation in both directions are merged, because EB
re-reflux can redistribute a cut-cell correction over a 3-by-3 neighborhood.
The planner rejects any result that cannot remain strictly inside the root
physical boundary.

A valid patch set contains separated, ratio-aligned children over one root.
Its composite integral is

`I = sum(uncovered root fluid volumes * U_c)`
`  + sum(all child fluid volumes * U_f)`,

so a root cell covered by any child contributes only through fine children.
Topology replacement applies

`average-down(old set) -> PCM(new set) -> retain(all old/new intersections)`.

The exact-overlap step uses global fine indices and may copy from any old child
to any new child, so movement, resizing, splitting, and repartition preserve
same-resolution data wherever physical fine cells coincide. Removal is the
empty new-set case. Candidate root and child arrays commit only after geometry,
finiteness, and active-cell EOS checks succeed.

For one coarse hydro interval, the root advances once. Every child then takes
`r` steps of `dt/r` with exterior states interpolated between the old and new
root states. Child `p` accumulates the exact coarse and fine face fluxes in its
own register `F_p`. The private root is corrected by each `F_p` in deterministic
child order and the complete set is average-downed after all corrections. The
required two-cell patch separation prevents redistribution neighborhoods from
coupling two independently refluxed interfaces.

With chemistry enabled, the patch-set Strang interval is

`R_root,children(dt/2) -> H_patch-set(dt)`
`  -> R_root,children(dt/2) -> average-down(all children)`.

Every reaction call masks covered EB cells. Coarse state, all fine states, all
temperatures, and patch metadata are private candidates until every reaction,
hydro, reflux, EOS, and synchronization stage succeeds. This kernel currently
has no molecular-transport stage.

## Public reactive EB multipatch lifecycle

Multipatch mode begins from the configured seed rectangle so initialization
has a valid hierarchy even when initial tagging is disabled. With initial
tagging enabled, and later at every `regrid_interval`, the root temperature
produces a new collection plan. An identical collection is a no-op. Empty tags
retain the current set unless `remove_fine_patch_when_untagged` is enabled; in
that case the empty-set topology transaction conservatively returns all child
data to the root.

The public coarse timestep is

`dt = min(dt_root, min_p(r_p * dt_child,p), remaining_time)`.

The accepted-step transaction advances that interval with the set-wide Strang
operator. Time, minimum-timestep, step, and regrid counters change only after
the physics or topology operation commits. The composite initial and final
integrals use every uncovered root cell and every active child exactly once.

The executable writes the synchronized root through the existing EB CSV path.
For child index `p`, it inserts `_patchNNNN` before the configured fine-output
extension and writes that child's actual geometry and state. Child order is
the deterministic collection order. Multipatch checkpoint/restart uses a
distinct schema because single-patch checkpoint schema one cannot represent a
child collection.

The patch-set schema records a strict signature for the mesh, EB geometry,
refinement, chemistry, hydro, redistribution, collection planner, and regrid
cadence. It then records time, counters, minimum accepted timestep, base
density, the root payload, and each ordered child's actual bounds and payload.
Final time, maximum steps, output paths, and checkpoint scheduling remain
restart-adjustable. A reader rebuilds every geometry, treats conserved state as
authoritative, recovers active temperatures through the EOS, validates the
complete separated set and terminal marker, and commits no output until all
children succeed.

Scheduled patch-set checkpoints are written only after an accepted Strang
interval and any due periodic regrid commit. A stop-after-write exits with the
stored time. Restart resumes from the stored step and regrid counters, so the
next periodic topology evaluation has the same cadence as an uninterrupted
run.

## Configured outflow-side EB patch

A configured single fine rectangle may include the first or last root cell in
either coordinate direction when that domain side uses outflow conditions.
For a fine face on a coarse/fine interface, the exterior conserved state and
temperature continue to come from the adjacent root cell, interpolated between
the root states at the beginning and end of the coarse interval. For a fine
face coincident with a physical side, there is no adjacent root cell. The
exterior state is therefore copied from the current fine boundary cell, giving
the same zero-normal-gradient closure as the qualified outflow boundary.

The fine state used for this physical-side closure advances after each
substep, rather than being frozen at the start of the coarse interval. Exterior
construction validates its dimensions, finiteness, and positive temperature
before advancing. Because a physical side has no uncovered coarse neighbor,
the EB flux register skips it and refluxes only the remaining coarse/fine
interfaces.

From `0.97.0`, the temperature tagger evaluates every active root cell. An
interior cell compares its four cardinal neighbors as before; a boundary cell
uses only cardinal neighbors that exist inside the domain, and all covered
neighbors remain excluded. Plan bounding, buffering, and minimum-size growth
use the complete root index range. Multipatch flood fill likewise traverses
boundary tags, preserves deterministic scan order, and coalesces candidates
under the same two-cell redistribution separation rule. Both single-patch and
patch-set topology transactions can therefore create, move, or resize an
outflow-side fine rectangle. Patch-set hydro uses the evolving child state on
each physical side and root-time interpolation on every remaining interface.

## Static three-level EB synchronization

For two strictly nested ratio-aligned rectangles, the three-level composite
integral is

`I = sum_root\middle(V U) + sum_middle\finest(V U) + sum_finest(V U)`.

Synchronization is ordered deepest first. The finest field is
volume-weighted onto a private middle candidate with the qualified two-level
EB average-down, and that candidate is then restricted onto a private root
candidate. Covered parent volume receives no contribution, while uncovered
parent cells remain bitwise unchanged. The reactive form repeats the same
ordering with conserved state and EOS temperature recovery at each parent.
Input finiteness, geometry, dimensions, and positive temperatures are checked
before work; failure at either stage returns the original root and middle
fields. This milestone does not advance, subcycle, or reflux the three-level
hierarchy.

## Three-level reactive EB hydrodynamics

Let the root-to-middle and middle-to-finest refinement ratios be `r1` and
`r2`. One accepted root interval applies

`root(dt) -> r1 * [middle(dt/r1) -> r2 * finest(dt/(r1*r2))]`.

The root prediction supplies time-interpolated exterior state to every middle
substep. Each middle prediction similarly supplies the finest exterior during
its nested substeps. PCM samples the parent at the start of a substep;
characteristic PLM samples its midpoint, matching the qualified two-level
temporal closure.

An independent flux register accumulates root/middle flux mismatch over the
whole root interval. A fresh middle/finest register closes each middle
interval: it refluxes first, then volume-weighted average-down synchronizes the
middle before the next substep. After all middle work, the outer register
refluxes root and middle, and the three-level EOS transaction synchronizes
finest to middle to root. Every field is private until all stages succeed.

The finest rectangle must be separated from the middle boundary by two cells,
and every face on its coarse/fine interface must have unit open-area fraction.
This regular-interface condition avoids silently applying the two-level EB
register to an EB-cut nested interface, which needs a dedicated geometric
correction. Root and middle geometry may still contain regular, cut, and
covered cells.

## EB-cut nested-interface conservation closure

For an EB-cut finest interface, the inner two-level update first performs its
normal flux-register reflux and average-down. Let `I0` be the middle/finest
composite integral before the middle interval, and let `B` be the signed,
time-integrated Cartesian flux through the exterior boundary of the middle
mesh. The conservative target for components without an EB wall source is

`Itarget = I0 + B`.

The closure compares this target with the post-sync composite integral.
Density, total energy, and every species residual are retained; momentum is
excluded because the embedded slip wall supplies a physical pressure force.
The sum of species residuals must agree with the density residual to scaled
roundoff. A final species component absorbs only that roundoff difference so
each corrected recipient retains local density/species closure.

The residual conserved density is divided by the total fluid volume of active
middle cells outside the finest rectangle and added uniformly to those cells.
This recipient set keeps the correction authoritative after finest-to-middle
average-down and away from the nested interface by the established two-cell
separation. Every corrected cell undergoes EOS temperature recovery. Invalid
composition, energy, or temperature rejects the complete hierarchy. This is a
qualified global multilevel conservation closure, not parity with PeleC's
locally resolved multilevel EB redistribution stencil.

## Three-level reactive EB Strang composition

For a root interval `dt`, the static three-level reactive driver applies

`R(dt/2) -> H3(dt) -> R(dt/2) -> A3`,

where `R` advances chemistry independently in every active root, middle, and
finest cell, `H3` is the recursively subcycled three-level EB hydro
transaction, and `A3` is reactive deepest-first average-down. Covered cells
are excluded from each reactor call.

The final `A3` operation is required even though hydro already ends
synchronized: the second reaction half-step also advances geometrically
overlapped parent cells, which are not part of the composite solution. The
finest result therefore replaces its middle footprint before the resulting
middle field replaces the root footprint, with EOS temperatures recovered at
both parents. All work uses private candidates; any rejected reactor cell,
hydro stage, nested EB conservation correction, or final EOS recovery leaves
every input state and temperature unchanged.

## Static three-level EB AMR time loop

Let `dt0`, `dt1`, and `dt2` be the active-cell stability limits computed on
the root, middle, and finest meshes. With refinement ratios `r1` and `r2`, the
root interval is

`dt = min(dt0, r1*dt1, r1*r2*dt2, final_time - time)`.

The ratio factors convert each child's substep limit to the equivalent root
interval. A successful interval advances through the complete three-level
Strang transaction before time, step count, and minimum accepted timestep are
updated. Initialization uses PCM prolongation first from root to middle and
then from middle to finest, so all levels begin EOS-consistent and synchronized.

Namelist validation always rejects multipatch ownership and colliding output
paths when `three_level_enabled` is selected. Static and dynamic three-level
modes use distinct checkpoint magics so their topology contracts cannot be
misinterpreted.

## Static three-level EB AMR checkpoint transaction

The static hierarchy uses a dedicated formatted stream whose magic and schema
are distinct from the single-patch and patch-set formats. It stores ordered
species names, root geometry and EB parameters, hydro and chemistry controls,
both nested rectangles and refinement ratios, accepted time, minimum timestep,
step count, base density, and all root, middle, and finest conserved and
temperature fields. Final time, maximum steps, output paths, and checkpoint
scheduling remain restart-adjustable.

Restart reconstructs all three meshes from the validated input topology,
reads every payload into private candidates, and rejects any schema, mechanism,
physics, topology, dimension, finiteness, or terminal-marker mismatch. Stored
temperature is checked for finiteness but is not trusted as thermodynamic
authority: active-cell temperature is recovered from conserved state through
the EOS on every level. Only a completely valid hierarchy is published.
Periodic and final checkpoints are written after an accepted Strang interval;
stop-after-write exits only after the committed stream is complete.

## Tag-driven three-level finest regridding

Dynamic three-level mode retains the configured root-to-middle patch and plans
only the middle-to-finest rectangle. Temperature-gradient tags are evaluated
on active middle cells. The outer two-cell band is excluded before buffering
and minimum-size growth, so any accepted plan satisfies the multilevel
redistribution stencil margin. If no interior tags remain, the existing
finest patch stays active.

For a changed rectangle, the old finest state is first volume-weighted onto a
private middle candidate. The replacement EB geometry is built from the new
middle-index bounds and PCM-prolonged from that synchronized candidate. Cells
overlapping the old and new finest rectangles retain their prior fine state
and temperature after geometry consistency checks; newly refined cells keep
the prolonged values. EOS recovery validates all active replacement cells
before the middle state, finest fields, geometry, and patch metadata are
published together.

The same transaction may run after initialization and after any accepted root
step selected by `regrid_interval`. The initial composite integral is measured
after an initialization-time topology change. Finest removal, changing the
root-to-middle patch, and sibling finest patches remain outside this mode.

## Dynamic three-level checkpoint transaction

The dynamic format extends the three-level compatibility contract without
changing the static stream. It stores the actual finest rectangle after any
committed regrid, the accepted regrid count, and every control that can change
future tag planning: interval, initialization flag, relative and absolute
temperature thresholds, scale floor, buffer, minimum dimensions, and maximum
patch gap.

Restart first validates the fixed root-to-middle rectangle and all physics and
regrid controls. It then validates the stored finest bounds against the
two-cell middle margin, builds that geometry rather than the namelist seed,
loads all three fields into private candidates, and recovers active
temperatures from conserved state. Time, minimum timestep, step count, regrid
count, topology, and fields publish together only after a complete end marker.
Checkpoint scheduling and output paths remain restart-adjustable.

## Single-level EB molecular transport

The regular mixture transport operator supplies Cartesian face-centered
viscous stress, Fourier heat flux, mixture-averaged species diffusion,
barodiffusion, correction velocity, and species enthalpy flux. These fluxes
are interpolated to the EB face centroids and multiplied by the nodal-geometry
open-face fractions. No diffusive source is added on the cut face itself, so
the qualified embedded wall is adiabatic, free slip, and impermeable to every
species.

For an active cell with fluid volume `V`, the conservative transport source is
the negative open-area flux sum divided by `V`. Before divergence, each
species' outgoing integrated flux is compared with 90 percent of its
fluid-volume inventory. The minimum adjacent-cell factor scales the complete
coupled face flux, retaining species/enthalpy and viscous-work consistency.

Each Euler transport stage passes that source through the established EB
StateRedist transaction and recovers active temperatures from conserved state.
Two such stages form SSPRK2. The public reactive EB step uses
`R(dt/2) -> T(dt/2) -> H(dt) -> T(dt/2) -> R(dt/2)` and publishes only the
fully valid result. The timestep is the minimum of the hyperbolic and mixture
transport stability limits.

## Two-level EB AMR molecular transport

For a coarse transport Euler interval `dt`, the coarse level advances once and
its EB face-centroid diffusive fluxes enter the interface register with weight
`dt`. The fine patch advances `r` Euler substeps of size `dt/r`. Before each
substep, its exterior conserved state is interpolated in time between the
coarse interval endpoints and converted to primitive transport data. Physical
outflow sides retain zero-gradient fine data. Each fine flux enters the same
register with weight `dt/r`.

After the fine substeps, EB-aware reactive reflux applies the integrated
coarse/fine flux mismatch and StateRedist handles any cut-cell correction.
Reactive average-down then makes the fine solution authoritative under the
patch. A second synchronized Euler transaction followed by the arithmetic RK
average forms SSPRK2. Temperature is recovered through the mixture EOS before
the final average-down. Consequently each Euler endpoint and the completed RK
step preserve the composite conserved quantities.

The root transport limit is

```text
dt_root <= min(dt_transport,coarse, r dt_transport,fine),
```

because the fine operator takes `r` temporal substeps during each root
interval. The public two-level driver inserts half intervals of this operator
around hydrodynamics and inside the chemistry half steps. Missing or mismatched
transport data, invalid exterior state, reflux failure, or EOS failure rejects
the whole hierarchy transaction. Three-level and sibling-patch transport are
not part of this milestone.

## Arbitrary-depth reactive EB patch-tree CFL reduction

For patch `p` on runtime level `l`, the established active-cell EB kernel
computes

```text
dt(l,p) = CFL / max_active[(|u| + c)/dx + (|v| + c)/dy].
```

If relation `q` has refinement ratio `r(q)`, define the number of subcycles from
level `l` to the root as

```text
R(1) = 1,  R(l) = product(q=1..l-1, r(q)).
```

The hydro-only stable root interval is therefore

```text
dt_root = min_over_all_nodes[R(l) dt(l,p)].
```

The traversal includes every branch and ignores covered cells through the
shared single-node kernel. A fully covered node contributes no constraint; an
entirely inactive tree, nonfinite control, or failed active-node conversion
rejects with deterministic zero output. Selection does not mutate topology,
state, or temperature.

For full physics, the same traversal also evaluates the explicit viscosity,
conduction, and species-diffusion limit on each active node. The root interval
is the minimum of the scaled hyperbolic and active transport limits over all
nodes. Barodiffusion changes the species flux but adds no independent
diffusivity constraint.

## Arbitrary-depth reactive EB patch-tree hydrodynamics

Let a parent node advance over interval `dt` and let its child relation have
ratio `r`. The parent is advanced once with the qualified EB level operator.
Every child is then advanced recursively `r` times with interval `dt/r`.
Exterior state for child substep `s` is interpolated between the parent's
start and uncorrected end states at

```text
alpha = (s - 1)/r                 for PCM,
alpha = (s - 1/2)/r               for characteristic PLM.
```

One register per child accumulates the parent flux with weight `dt` and each
child flux with weight `dt/r`. After all substeps, deterministic child-order
reflux and reactive average-down make the finer state authoritative. The same
operation recurses without a compile-time level limit or single-parent chain
assumption.

For every node with children, subtree conservation is checked against flux
through that node's outer EB boundary. Density, total-energy, and species
residuals are closed only on active parent cells not represented by a child;
species correction is constrained to equal the density correction, and all
corrected temperatures are recovered through the NASA7 EOS. A final
deepest-first average-down preserves the composite state after ancestor reflux.
The whole hierarchy commits only after every recursive operation succeeds.

## Arbitrary-depth reactive EB patch-tree chemistry

For each tree node `(l,p)`, construct the logical chemistry mask

```text
active(i,j) = cell_type(i,j) /= covered.
```

The established constant-volume 2D reaction integrator advances only those
active cells. Chemistry is cell-local, so every patch advances over the same
physical reaction interval independent of level and refinement ratio. After a
standalone chemistry transaction, deepest-first average-down restores the
composite representation before the candidate commits.

The coupled tree operation uses Strang ordering

```text
chemistry(dt/2) on every patch
recursive EB hydrodynamics(dt)
chemistry(dt/2) on every patch
deepest-first synchronization
```

All stages mutate one private candidate. Consequently, a late chemistry or
hydrodynamics rejection cannot expose a partially reacted hierarchy, and
per-level chemistry and hydro call counts are published only with the final
valid state. Molecular transport is not part of this `R-H-R` operation.

## Arbitrary-depth reactive EB patch-tree molecular transport

For one Euler stage, a node over interval `dt` evaluates the established
mixture diffusive flux and EB redistribution operator once. A child relation
with ratio `r` then advances recursively for `r` intervals of `dt/r`, using

```text
alpha = (s - 1)/r
```

to interpolate its exterior between the parent start and Euler-end fields.
One diffusive EB flux register per child accumulates the parent flux with
weight `dt` and each child flux with weight `dt/r`. Reflux, average-down, and
outer-flux subtree closure restore the composite conserved representation.

The public second-order update is SSPRK2 over the complete tree:

```text
U(1)     = Euler_tree(U(n), dt)
U(2)     = Euler_tree(U(1), dt)
U(n + 1) = 0.5 [U(n) + U(2)].
```

The blend is applied to every runtime node, followed by active-cell EOS
temperature recovery and deepest-first synchronization. Both Euler trees and
the blend remain private until the final candidate validates; the minimum
positivity-limiter theta and actual recursive node-call counts follow the same
commit boundary.

## Arbitrary-depth reactive EB patch-tree full physics

The complete split interval is

```text
R(dt/2) -> T_SSPRK2(dt/2) -> H(dt) -> T_SSPRK2(dt/2) -> R(dt/2).
```

Here `R` advances every active runtime patch over the same reaction interval,
`T_SSPRK2` performs two recursively subcycled diffusive Euler trees, and `H`
performs one recursively subcycled hyperbolic tree. Each operator retains its
own deepest-first synchronization and conservation closure, while all
intermediate hierarchies remain inside one outer candidate.

The final chemistry stage is followed by deepest-first synchronization and
tree validation. Only then are conserved state, temperature, the minimum theta
from both transport half-steps, and the chemistry/transport/hydro node-count
vectors committed. This outer transaction makes a rejection after any valid
prefix observationally equivalent to no attempted step.

## Arbitrary-depth reactive EB patch-tree time loop

At accepted time `t`, the public clock recomputes the combined all-node stable
interval and clips it to the requested target:

```text
dt = min(dt_hydro_and_transport(tree), final_time - t).
```

One complete `R-T-H-T-R` operation runs on a private tree candidate. Only
after that candidate succeeds are the tree, `t + dt`, total step count,
minimum accepted interval, minimum transport theta, and per-level physics
counts published. The next stability calculation therefore always observes
the last accepted state.

If the first step rejects, state and all clock/accounting outputs remain at
their input values or deterministic neutral outputs. If a later step rejects
or the maximum step count is reached, earlier accepted steps remain visible
and the failed candidate does not. Reaching the target within the floating-
point time tolerance assigns the requested final time exactly.

## MPI arbitrary-depth EB patch-tree ownership

For runtime node `(l,p)` with allocated cell count `N(l,p)`, the initial
distribution assigns work

```text
W(l,p) = N(l,p) R(l)^e,
```

where `R(l)` is the cumulative refinement product and the caller selects
`e = 0`, `1`, or `2`. Nodes are visited in deterministic level/patch order and
assigned to the rank with the smallest accumulated work, breaking ties toward
the lower rank. The distribution stores exact per-rank cell, node, and work
totals and validates them against every owner entry.

The first publication operation keeps a complete candidate tree on every rank.
Only the assigned owner supplies each node's conserved state and temperature;
two ordered broadcasts populate that node in the candidate. Collective input
and final-candidate validation surround the complete traversal, so no partial
publication is committed. This operation establishes numerical ownership but
does not yet remove nonowner field allocation or route recursive physics.

## MPI sparse arbitrary-depth EB patch-tree storage and migration

The sparse tree separates replicated metadata from numerical fields. Every
rank stores the same topology and owner map, but node `(l,p)` allocates state
and temperature exactly when the local rank equals its owner. Validation
checks both directions of this invariant, along with exact geometry-derived
array shapes, finite fields, and positive temperature.

Materialization is an explicit compatibility boundary: allocate a private
complete tree, copy each owner's node into it, broadcast that node's state and
temperature, validate, then publish. Repartitioning instead allocates a private
sparse tree from the new owner map. Unchanged owners copy locally; a changed
owner pair performs direct point-to-point state and temperature transfer. The
accepted sparse tree remains untouched until every candidate node validates
against the new distribution.

## MPI owner-local arbitrary-depth EB patch-tree timestep

For every active node owned by a rank, compute the same hyperbolic and optional
explicit-transport limits as the serial patch tree. If `R(l)` is the cumulative
refinement product from the root to level `l`, each local node contributes

```text
dt_root(l,p) = R(l) min(dt_hydro(l,p), dt_transport(l,p)).
```

The rank-local minimum is reduced with `MPI_MIN`. A rank that owns no active
node contributes `huge`, so it cannot constrain a valid result; a separate
active-node sum rejects an entirely covered tree. Collective preflight requires
identical CFL values, transport-enable flags, and species count, plus valid
owner-only fields and a replicated distribution descriptor. Thus timestep
selection neither broadcasts nor materializes numerical fields.

## MPI owner-local arbitrary-depth EB patch-tree chemistry

Each sparse owner applies the established cell-local constant-volume reactor to
its node with the EB active mask, then reconstructs active temperatures from
the conserved state. All ranks visit nodes in the same deterministic order and
collectively accept after each owner operation. Chemistry counts therefore
publish only for the final committed sparse candidate.

After all nodes react, relation `L` is synchronized from deepest to shallowest.
For child `c` and parent `p`, the child state is already authoritative after
all deeper descendants have restricted into it. If `owner(c) != owner(p)`, the
child state is sent directly once to `owner(p)`; otherwise no communication is
needed. The parent owner applies the serial EB average-down operation in child
order and retains the updated parent state and temperature. No child
temperature or complete tree is transmitted.

## MPI sparse arbitrary-depth EB composite integrals

For a selected node, construct a Boolean mask from all of its direct child
rectangles. The node owner contributes

```text
I_local = sum(unrefined cells) alpha(i,j) U(i,j) dx dy,
```

where `alpha` is EB volume fraction. Every rank then follows every descendant
relation, but only the owner of each visited node evaluates cells. A single
communicator `MPI_SUM` produces the composite conserved vector; a companion
integer sum verifies that at least one owner node contributed and supplies
diagnostic accounting.

The complete-tree operation selects the root. A subtree selection is accepted
only when level, patch, and vector extent agree on every rank. This makes the
reduction suitable for recursive conservation closure without introducing a
replicated numerical-tree boundary.

## MPI owner-local arbitrary-depth EB patch-tree hydro

For each node invocation, only `owner(level,patch)` advances the conserved
state and temperature and retains the node's Cartesian EB fluxes. For every
child edge, the parent start/end boundary context is transferred once when
ownership differs. The child then performs `r` recursive substeps using the
same time interpolation as the serial tree. After each substep its fine fluxes
move directly to the parent owner for accumulation into that edge's register.

After all child substeps, each register is consumed in child order. Across an
owner boundary, current child state moves child-to-parent, the parent applies
reflux to its own current state, and corrected child state returns. One later
child-to-parent transfer supplies ordered average-down. Thus a remote edge
executed by one parent invocation produces

```text
r + 4
```

grouped direct transfers: one context, `r` fine-flux payloads, two reflux
payloads, and one average-down payload. Shared-owner edges produce none.

Every refined node uses owner-local subtree reductions for its before/after
composite integral and applies the same boundary-flux conservation residual to
unrefined active parent cells. The public sparse operation commits only after
the complete root recursion and candidate validation succeed.

## MPI owner-local arbitrary-depth EB patch-tree transport

Each recursive transport Euler stage uses the hydro ownership route with the
transport flux/RHS and StateRedist kernels. A remote parent/child edge sends
one time-interpolated exterior context, returns one fine diffusive-flux payload
for each of the `r` child substeps, performs the two-message reflux round trip,
and sends corrected child state once for average-down. Its grouped transfer
count is therefore again `r + 4` per parent invocation.

Let `U^0` be the accepted sparse tree and let `E(U)` denote one owner-local
recursive Euler stage. The SSPRK2 transaction computes

```text
U^1 = E(U^0),
U^2 = E(U^1),
U^(n+1) = 1/2 (U^0 + U^2).
```

The final blend and EOS temperature recovery occur independently on every node
owner. A deepest-first direct average-down follows the blend, adding one remote
transfer for every distinct-owner relation. The public limiter is the
communicator minimum over both stages. All fields and diagnostics commit only
after final sparse validation.

## MPI owner-local arbitrary-depth EB full physics

For accepted sparse tree `U^n`, timestep `dt`, reaction operator `R`, SSPRK2
transport operator `T`, and recursively subcycled hydro operator `H`, compute

```text
U* = R(dt/2) T(dt/2) H(dt) T(dt/2) R(dt/2) U^n.
```

Disabled chemistry omits both `R` applications. Disabled explicit transport
makes each `T` application an accepted no-op with unit limiter. Every active
operator consumes and returns only owner-local sparse fields through its
qualified communication schedule.

The outer operation accumulates chemistry-node, transport-Euler, and hydro-node
advances independently, along with their three transfer categories. Its
minimum transport limiter is the minimum from both half-steps. The accepted
tree is assigned from `U*` only after every stage and final sparse validation
succeed; failures expose neither a partial split state nor partial accounting.
