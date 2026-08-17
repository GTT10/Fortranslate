# Numerical methods

## Governing equations

The current solver advances

\[
\partial_t U + \partial_x F(U)=0,
\qquad
U=(\rho,\rho u,\rho v,\rho w,\rho E)^T,
\]

with a constant-`gamma` ideal-gas closure.

## PeleC limited slopes

For cell values `q_{i-2},...,q_{i+2}`, define the centered differences and monotonized bounds used in PeleC. With order 2, the neighboring limited slopes are zero and the returned slope reduces to the limited centered expression.

With order 4, limited neighboring slopes `d_{i-1}` and `d_{i+1}` are first constructed, then

\[
d_{\mathrm{temp}}=
\frac{4}{3}d_{\mathrm{cen}}
-\frac{1}{6}\left(d_{i-1}+d_{i+1}\right),
\]

and

\[
d_i=f_i\,\operatorname{sign}(d_{\mathrm{cen}})
\min\left(d_{\mathrm{lim}},|d_{\mathrm{temp}}|\right).
\]

Here `f_i` is the flattening coefficient. `plm_order=2` and `plm_order=4` are both retained for differential testing.

## Shock flattening

The implemented flattening coefficient follows the regular-cell one-dimensional PeleC logic. It uses pressure values through `i+3`, normal velocity through `i+1`, and the constants

```text
shock threshold = 0.33
zcut1           = 0.75
zcut2           = 0.85
```

A pressure jump is flattened only when it is sufficiently strong and the local velocity field is compressive. A second, pressure-jump-directed shifted test is combined with the local test. The result is

[[
f_i=1-\max(\chi z,\chi_2 z_2),
\qquad 0\le f_i\le 1.
\]

Smooth regions return one. A strong compressive step may return zero. Expansion waves are not flattened.

## Characteristic tracing

Primitive slopes are projected into left acoustic, contact, two shear, and right acoustic waves. The waves are traced with speeds `u-c`, `u`, and `u+c` over the current `dt/dx`, and then transformed back into time-centered interface states. Invalid reconstructed states fall back to the corresponding cell-centered state.

## Sedov-type strong-blast case

The regression domain is `[0,1]` with a blast centered at `0.5`. Density is initially uniform. Pressure is `100` inside radius `0.0125` and `1e-5` outside. The simulation advances to `t=0.02` on 800 cells with fourth-order characteristic PLM and flattening.

This is a symmetric planar blast, not the spherical analytic Sedov-Taylor similarity problem. It is used as a strong-shock verification gate for:

- finite and positive states;
- reflection symmetry of density and pressure;
- antisymmetry of velocity;
- mass, momentum, and energy conservation;
- shock-radius and field-signature stability.

## Verified results

The fourth-order periodic entropy-wave test gives density L1 errors

```text
nx=40    3.5852193457e-4
nx=80    7.9925449806e-5
nx=160   1.7613156678e-5
```

with observed orders `2.1653` and `2.1820`. The overall scheme remains second order in time, as expected.

The 800-cell blast completes in 1052 steps. Representative final metrics are:

```text
minimum density       1.4352255263e-1
maximum density       5.0001945148
minimum pressure      1.0e-5
maximum pressure      1.3654133597e1
shock radius          1.31875e-1
mass error            1.58e-14
energy error           3.29e-14
```

## Scope limitations

The implementation does not yet include fourth-order multidimensional stencils, embedded-boundary slope fallbacks, general-EOS internal-energy characteristic amplitudes, species/passive-scalar tracing, rotating frames, or transverse Godunov corrections.
