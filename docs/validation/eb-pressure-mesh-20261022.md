# Native MOL thermal-contact mesh diagnostic, 2026-10-22

This batch follows the completed recursive pressure integration and its
50/25 ns timestep comparison. It tests the remaining nonisothermal pressure
error at 8, 16 and 32 cells per axis on each of the two AMR patches. It is a
new hydro-only diagnostic, not a rerun of the reacting PeleC qualification.

## Matched inputs and operator

`tools/check_eb_pressure_mesh.py` reads the saved 8/16/32 imported states.
All three source reports identify the same profile SHA256:
`cb9d6f058f23430399d8d75b2448ace3897b4d716b62eb5ca9b4dbbc556334ef`.
The domain is 1 cm, the planar embedded boundary remains at z=0.0033 m,
and the refinement ratio is two. The cut-cell fraction changes with grid
spacing; this is fixed physical geometry, not a fixed cut fraction.

Density, species densities, held energy and temperature are preserved from
each saved import. Momentum is replaced with rho*(10,0,0), and total energy
with held energy plus kinetic energy. The imported reaction source is unused
by the hydro-only audit. Initial EOS pressure must be within 1e-7 Pa of
101325 Pa. The saved imports and the canonical reference inputs are untouched.

The same frozen production-core probe used by the recursive integration
record runs four 50 ns root steps, two MOL iterations, recursive clocks,
no chemistry and no molecular transport. Both policies use consistent velocity
and separate momentum; only paired pressure differs. The comparison endpoint
is always 200 ns. Composite budgets include returned numerical-energy
corrections and must remain within 1e-12. Full inventory CSVs remain available.

## Execution and interpretation

Artifacts and atomic progress: `build/eb-pressure-mesh-20261022/record.json`.
Driver log: `build/eb-pressure-mesh-20261022.driver.log`.
The record identifies the driver and active child PID, exact executable,
commands, inputs, hashes, durations and completed results. A frozen source
manifest, runner and input lineage record are kept beside it.

The batch **completed** all three resolutions, both policies and all four
steps per trajectory. There were no failures, skipped grids or timeouts.
Elapsed times were 41.6, 121.2 and 521.3 seconds (about 11.4 minutes total).
The 16-cubed pilot forecast 1454.9 seconds for 32 cubed including margin;
the actual final grid completed within the 3600-second limit.

## Results at 200 ns

| Cells per patch axis | Velocity-only max pressure error (Pa) | Paired max pressure error (Pa) | Paired maximum normalized budget residual |
| --- | ---: | ---: | ---: |
| 8 | 6.2389464e-2 | 1.2712569e-3 | 1.03161e-14 |
| 16 | 1.8680804e-2 | 1.3900369e-4 | 5.34374e-14 |
| 32 | 4.7381278e-3 | 1.4555524e-5 | 3.10308e-13 |

Observed orders log2(error_coarse/error_fine) are 1.740 and 1.979 for
velocity-only, and 3.193 and 3.255 for paired pressure. The paired error
falls about 87.3-fold from 8 to 32; at 32 it is about 325.5-fold smaller
than velocity-only. Initial pressure error is at most 4.37e-11 Pa.
All normalized budgets pass the unchanged 1e-12 threshold, including the
32-cubed accumulation residual. Returned signed numerical-energy corrections
are included; the paired final corrections are 8.26e-18, 5.94e-18 and
2.93e-18 J. These are signed global corrections, not absolute local norms.

`audit.py` independently recomputes all 420 component-inventory residuals
from initial/final inventories and physical/numerical changes, checks all
six trajectory summaries, recomputes both policy orders, and verifies frozen
binary, runner and imported-input hashes. `audit.json` records the passed
checks and artifact hashes. Provenance, Python syntax, project contract and
`git diff --check` also passed. No solver code changed for this batch.

## Qualification boundary

This closes the bounded question of whether the residual pressure error
shrinks under spatial refinement for this nonisothermal moving contact:
it decreases on both refinements with closed correction-accounted budgets.
Together with the prior coarse-grid timestep halving, this supports spatial
error as the dominant residual in this control.

The measured 3.2 order is for one pressure-error observable at one short
endpoint, with fixed dt and a shared tabulated initial profile. It is not a
claim of third-order accuracy of the solver, a joint space-time asymptotic
study, or convergence of all flow variables. Cut fractions change with the
fixed physical plane. The saved source reports include a failed reacting
reference comparison; only their initial fields are reused here, and that
historical failure is not reclassified as passed.

Pinned reacting PeleC parity and broad physical validation remain open.
The paired policy intentionally changes the reference spatial algorithm.
No production defaults or frozen reference acceptance values were changed.
This batch has no remaining running job.
