# Completed reacting EB-AMR mesh comparison, 2026-10-21

The 8-, 16- and 32-cubed root/fine diagnostic trajectories all satisfy the
unchanged field, species, inventory and energy-budget limits at 0, 50 and
100 ns. These use explicit density chemistry tolerance `1e-18 kg/m3`,
relative tolerance `1e-12`, PeleF's opt-in converged caloric inverse, the
native activation-basis mechanism and the native held-EOS diagnostic.
The original pinned solver and canonical mechanism remain separate authorities.
This closes this configured three-grid trajectory comparison, not C2 as a whole.

## Completed results and cost

| Root/fine patch size | Native wall time, s | PeleF wall time, s | Maximum temperature difference, K | Maximum total-energy difference, J/m3 |
| --- | ---: | ---: | ---: | ---: |
| 8 cubed | 17.475 | 41.404 | 3.988e-5 | 1.202e-5 |
| 16 cubed | 102.508 | 224.198 | 1.182e-5 | 3.030e-5 |
| 32 cubed | 654.584 | 1540.658 | 2.732e-5 | 3.574e-4 |

The maximum energy-budget relative residual across all three grids/times is
`6.86e-17`. Relative H2-consumption differences are about `2.24e-8`.
The new 16/32 runs and conservative mesh analysis take 2569.38 s (42.82 min)
including exports/checking. The 8-cubed result was reused from the preceding
completed implementation check. These are Debug, one-thread measurements,
not Release throughput or parallel-scaling results.

Artifacts are in `build/eb-eos-mesh-study-20261021/`: `record.json` lists all
three exact reports, solver/checker commands and timings; `mesh-comparison.json`
records fluid-volume restriction residuals, initial/endpoint/increment norms
and observed orders. `job/result.json` records wrapper exit 0, and all solver,
export and comparison exits are zero. The checker and analyzer used at run
time are preserved as `*_at_run.py` before subsequent input-tool changes.

## What the observed orders mean

At a common time step, the observed orders of the evolution increments
(`state(100 ns)-state(0)`) are:

| Quantity | PeleF | Native held EOS |
| --- | ---: | ---: |
| Density | 1.356 | 1.356 |
| x-momentum | 1.518 | 1.518 |
| z-momentum | 1.094 | 1.094 |
| Total energy | 1.395 | 1.395 |
| Temperature | 1.444 | 1.444 |

These are observed differences, not proof of asymptotic second-order accuracy.
Endpoint temperature order is only 1.012, essentially the same as its
initial-data order. The fixed 33-point input table explains that initial
behavior: native TGReact linearly interpolates its tabulated values at each
Cartesian cell center. The chosen grid family changes its alignment with
the table's slope discontinuities. A saved-state audit reproduces every
initial temperature to `4.55e-13 K` and reproduces the measured initial order
1.01211179846 without advancing either solver.

## Separate smooth sampling control

`tools/prepare_tgreact_profile.py` verifies the frozen table against the
already documented formulas `T=2000+50*cos(2*pi*x/L)` and
`Y_H2=0.06+0.002*cos(2*pi*x/L)`, with `Y_O2=0.5`, balance N2 and zero trace
species. It writes a separate 129/257/513-point table and provenance manifest;
it refuses source overwrite, including a hard-link alias, and altered source
values. No original case input or acceptance hash is replaced.

With 129 points, every root/fine cell center on the 8/16/32 grids samples the
analytic values exactly. The independent initial-temperature restriction
order becomes 1.99653335489. This establishes the input-sampling effect only;
it does not establish time-evolution accuracy. The saved-field audit,
generated input/manifest and refusal checks are in
`build/eb-initial-mesh-audit-20261021/`. Python syntax, the two mesh-comparison
tests, `git diff --check` and the project contract pass after integration.

The diagnostic checker accepts `profile_points: 129` in mesh-study mode,
overrides only its working profile and native point count, and records the
sampling convention and actual profile hash. The mesh analyzer rejects mixed
sampling conventions or profile bytes. Its frozen-reference path is unchanged.

The separate 129-point 8/16/32 trajectory study completed with the same
time/chemistry/EOS controls. All three reports pass with zero violations;
the wrapper exits zero. Artifacts are in
`build/eb-smooth-mesh-study-20261021/`, including exact source identities,
commands and report paths in `record.json`, and orders in `mesh-comparison.json`.
Total elapsed time is 2622.05 s (43.70 min). Native/PeleF solver times are
17.720/41.645, 103.226/225.432 and 662.611/1522.184 s respectively.
Maximum temperature differences are 3.988e-5, 1.256e-5 and 2.590e-5 K;
the maximum relative energy-budget residual is 1.28e-16.

| Quantity | Endpoint order, PeleF | Evolution-increment order, PeleF |
| --- | ---: | ---: |
| Density | 1.980 | 1.356 |
| x-momentum | 2.106 | 1.518 |
| z-momentum | 1.093 | 1.093 |
| Total energy | 1.709 | 1.395 |
| Temperature | 1.992 | 1.444 |
| H2O | 1.990 | 1.990 |

The native held-EOS result has essentially the same increment orders.
Smoother initial sampling resolves the first-order initial-temperature effect,
but does not establish uniformly second-order evolution or isolate a
Fortran-specific defect.

An independent saved-field audit in
`build/eb-smooth-mesh-regions-20261021/{analyze.py,regions.json}` partitions
the fluid volume into six fixed physical regions and verifies that their
squared-error fractions sum to one. For the 16-to-32 evolution comparison,
the upper bulk holds 32.1% of density, 34.2% of energy and 36.3% of temperature
squared error. The lower refined slab holds 38.9% of x-momentum and 40.1%
of z-momentum error. The remaining discrepancy is therefore not confined to
the wall. These saved-data results do not separate limiter, refinement-interface,
boundary-compatibility and temporal effects; further asymptotic work remains.

Original-PeleC parity, asymptotic space/time convergence, full
configuration/install qualification and physical validation remain open.
