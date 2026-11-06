# Residual shear maxima under refinement — 2026-10-23

## Verdict

The connected-direction candidate's weak transverse Linf convergence is
localized at coarse/fine interfaces. The y-momentum and z-momentum maxima
have different geometric locations. A cut-cell-only correction cannot be
justified as a remedy for both: the paired z maximum is in a full-volume
fine cell at the intersection of x and z coarse/fine boundaries.
This is a saved-data diagnosis, not a new solver run or a proven causal fix.

## Evidence and reproducibility

`build/eb-shear-local-convergence-20261023/analyze.py` reads the completed
frozen candidate results in `build/eb-shear-tangent-convergence-20261022/runs/`.
It writes `analysis.json` with SHA256 identities and `locations.csv` with
physical coordinates, geometry flags, norms and growth measures. It checks
all saved norm/max-location pairs for equality of Linf and absolute signed
error, then examines 144 noninitial maxima: 3 grids × 4 policies × 4 times ×
3 variables (y momentum, z momentum and pressure). **All 144 are marked
coarse/fine adjacent.** The coordinates follow the probe's zero-origin nested
half-domain patch with refinement ratio 2, independently checked against
its initialization and spacing. Locations are Cartesian cell centers, not
cut-cell fluid centroids. Tied maxima retain the probe's first encountered cell.

No calculation was started, and no completion notification is needed for
this synchronous analysis. Python syntax, mechanism provenance, project
contract and `git diff --check` pass. This audit does not rerun the previous
conservation audit; the earlier 12-run maximum budget remains historical
candidate evidence (2.071e-13).

## Paired-policy endpoint maxima

The following are at 200 ns, in fine node 2. Momentum errors are conserved
momentum-density errors; the table does not report velocity errors.

| Root grid | y-momentum Linf | Cell (i,j,k) | Volume fraction | z-momentum Linf | Cell (i,j,k) | Volume fraction |
| --- | ---: | --- | ---: | ---: | --- | ---: |
| 8 | 1.958672e-3 | (1,4,6) | 0.72 | 9.422559e-4 | (1,5,8) | 1 |
| 16 | 2.035193e-3 | (1,8,11) | 0.44 | 9.270569e-4 | (1,9,16) | 1 |
| 32 | 1.842048e-3 | (1,16,22) | 0.88 | 7.929094e-4 | (1,17,32) | 1 |

All six maxima lie one half fine-cell width from the x coarse/fine interface
and one half fine-cell width from the shear extremum at y=0.0025 m.
The y maxima are cut cells at the embedded plane z=0.0033 m. Their changing
volume fractions (0.72, 0.44, 0.88) confound a pure fixed-stencil refinement
interpretation. The z maxima instead lie one half cell below the fine patch's
top at z=0.005 m, well away from the embedded plane. They have no domain-wall
adjacency flag. Thus the residual includes a coarse/fine intersection effect
in regular-volume cells; it is not confined to small-volume amplification.
The 32-grid y maximum also switches location between 50 and 100 ns, so maximum
histories are not necessarily a single-cell trajectory.

## Why L2 improves while Linf barely changes

For volume-weighted L2, define `f_eff=(L2/Linf)^2`. This is the normalized
integral of squared error divided by its maximum square. It is an equivalent
volume measure, not the measured support volume or a direct band-width estimate.

| Root grid | y-momentum f_eff | z-momentum f_eff |
| --- | ---: | ---: |
| 8 | 4.664936e-3 | 2.901715e-3 |
| 16 | 1.225734e-3 | 3.463993e-4 |
| 32 | 3.099838e-4 | 4.170212e-5 |

These shrink by roughly factors of 4 (y) and 8 (z) per refinement. This
supports concentration of the residual into progressively smaller weighted
regions even though its amplitude remains similar. It does not prove a
geometric codimension or establish pointwise convergence.

From 50 to 200 ns, `log(Emax(200)/Emax(50))/log(4)` is 1.999/2.003/1.925
for y momentum and 1.976/1.948/1.841 for z momentum. Pressure has exponents
0.974/0.940/0.795. The early behavior is consistent with pressure error
building first and subsequently driving transverse acceleration, as in the
previous stage diagnostics. These are growth exponents at fixed time step,
not temporal accuracy orders; they do not establish a pressure-gradient cause
at these particular maximum cells.

## Next numerical discriminator and qualification boundary

Before changing the limiter, capture directional face fluxes and pressure
reconstruction around the shear extremum at the x/z coarse-fine intersection
and at the x/cut-plane intersection. Compare the pressure differences on
opposite faces with transverse momentum-rate differences at matched stages.
The stored maxima alone contain no neighboring pressure or flux values, so
they cannot distinguish loss of smooth-extremum accuracy, halo interpolation,
or face synchronization as the specific remaining operator defect.

The production path matches the candidate at 16 cells/50 ns, as recorded in
[the integration comparison](eb-connected-shear-integration-20261023.md).
The 8/16/32 evidence here comes from the frozen candidate, not a new production
convergence campaign. The opt-in remains experimental; pointwise convergence,
general geometry and PeleC reference qualification stay open.
