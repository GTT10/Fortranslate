# Numerical methods

## Governing equations

The current solver advances the one-dimensional Euler system

\[
\partial_t U + \partial_x F(U)=0,
\qquad
U=(\rho,\rho u,\rho v,\rho w,\rho E)^T,
\]

closed by a constant-`gamma` ideal gas.

## Spatial update

For cell `i`,

\[
\frac{dU_i}{dt}=-\frac{F_{i+1/2}-F_{i-1/2}}{\Delta x}.
\]

The available Riemann fluxes are Rusanov and a single-species ideal-gas reduction of the star-state and wave-interpolation logic in PeleC `Source/Riemann.H`.

## Reconstruction paths

### `pcm`

Adjacent cell averages are sent directly to the selected Riemann solver.

### `plm`

Primitive variables

\[
q=(\rho,u,v,w,p)^T
\]

are reconstructed componentwise with minmod or MC slopes. SSPRK2 supplies second-order time integration.

### `pelec_plm`

The same limited primitive slope is decomposed into the five one-dimensional characteristic families: left acoustic, contact, two shear waves, and right acoustic. With local sound speed `c`,

\[
\alpha_- = \frac{\rho}{2c}\left(\frac{\Delta p}{\rho c}-\Delta u\right),
\qquad
\alpha_+ = \frac{\rho}{2c}\left(\frac{\Delta p}{\rho c}+\Delta u\right),
\]

\[
\alpha_0 = \Delta\rho-\frac{\Delta p}{c^2},
\qquad
\alpha_v=\Delta v,
\qquad
\alpha_w=\Delta w.
\]

The inverse mapping is verified independently:

\[
\Delta\rho=\alpha_-+\alpha_0+\alpha_+,
\]

\[
\Delta u=\frac{c}{\rho}(\alpha_+-\alpha_-),
\qquad
\Delta p=c^2(\alpha_-+\alpha_+).
\]

PeleC-style tracing then evolves these waves toward the two faces using speeds

\[
u-c,\qquad u,\qquad u+c
\]

and the current `dt/dx`. The resulting face states are passed to the selected Riemann solver. Because these states are already time-centered, the update is

\[
U_i^{n+1}=U_i^n-\frac{\Delta t}{\Delta x}
\left(F_{i+1/2}^{n+1/2}-F_{i-1/2}^{n+1/2}\right).
\]

## Positivity and fallback

Density and pressure slopes share a reduction factor if either extrapolated face approaches its configured floor. Failed primitive-to-conserved conversion falls back to the corresponding cell average. A timestep is rejected if the final interior state is non-physical.

## Boundary treatment

- `outflow`: copy the nearest interior state; suppress boundary-adjacent slopes.
- `periodic`: wrap both states and slopes.

## Verified accuracy

For periodic entropy-wave advection with `pelec_plm` and the PeleC-style Riemann solver:

| Cells | Density L1 |
|---:|---:|
| 40 | `3.4863582e-4` |
| 80 | `7.8837355e-5` |
| 160 | `1.7751480e-5` |

Observed orders are `2.1448` and `2.1509`.

For 400-cell Sod at `t=0.2`, density and pressure L1 errors are `1.2178e-3` and `7.1822e-4`.

## Scope limitation

The current characteristic implementation excludes PeleC fourth-order slope construction, flattening, general-EOS internal-energy characteristics, species/passive-scalar tracing, EB-specific paths, rotation terms, and multidimensional transverse corrections.
