# Nonuniform reacting EB-AMR reference comparison, 2026-10-19

The new comparison runs the production Fortran EB-AMR reaction/hydro/transport
integrator from the same initial gas state as unchanged pinned PeleC TGReact.
Both trajectories complete, but the field comparison **fails** the existing
regular-TGReact acceptance gates. This is a reproducible C2 discrepancy,
not a qualified EB-AMR trajectory result.

## Implemented comparison

`tools/reference_probes/eb_tree_reacting_probe.F90` is built as
`test_eb_tree_reacting_reference`. It imports all time-zero gas, temperature
and held internal-energy cells, independently initializes reaction history,
and advances the production recursive EB MOL driver. The required internal
energy flux storage is forwarded alongside UEINT. No production numerical
module was changed in this increment.

`tools/check_tgreact_eb_amr_reference.py` checks reference source hashes,
case inputs, physical times, hierarchy and analytic/exported volume fractions.
It imports Simple viscosity/conductivity/binary-diffusion fits from the exact
pinned LiDryer header, checking SHA256, Git blob and species order through the
existing importer. Thus the comparison uses the same molecular transport model.
Captured reference files can be reused only with matching source/input
contracts and file hashes. Native sources, golden data and tolerances were
not changed. The new probe is an explicit target, not another automatic CTest
matrix.

The [case](../../cases/tgreact_eb_amr_reference/README.md) has an 8³ root,
an 8³ fine patch over coarse indices 1:4 in all directions, a 1 cm cube,
plane z = 3.3 mm, periodic x/y and adiabatic slip z/EB boundaries. Initial
temperature span is 98.08 K and maximum speed 9.442 m/s. Three 100 ns root
steps use recursive MOL2, CVODE and molecular transport. The composite has
544 fluid cells: 352 coarse (48 cut) and 192 fine (64 cut). Source initialization
is independent; only time-zero gas/T/UEINT are required to match exactly.

## Final result and diagnostic control

Final corrected-boundary evidence is
`build/tgreact-eb-amr-reference-20261019/runs/matched_eb_amr_62p3fjsx/checks.json`.
The native run takes 5.78 s and the Fortran run 38.12 s. All four synchronized
outputs exist. The field gate has 21 violations across steps 1–3:

| Quantity | Maximum difference through 300 ns | Existing limit |
| --- | ---: | ---: |
| Temperature | 0.01498 K | 0.001 K |
| Density relative to field maximum | about 1.25e-4 | 1e-8 |
| x momentum density | 0.006622 kg/m²/s | 1e-6 kg/m²/s |
| Species mass fraction | 8.912e-8 | 1e-7 |
| H2 consumption relative difference | 2.242e-5 | 1e-3 |

The species and consumption gates pass. Native and Fortran consume respectively
`4.629508317e-5` and `4.629404527e-5` of initial H2. Independent maximum relative
mass/element drifts are `1.34e-15` native and `5.98e-14` Fortran; the Fortran
energy-budget residual after its recorded correction is below `9.41e-17`.
These conservation results do not cancel the spatial field failures.

The transport-off control (`matched_eb_amr_n7d04w93`) also completes and fails
with comparable density/tangential-momentum differences. Its first-step coarse
cut-cell x-momentum difference is `0.002466 kg/m²/s`, versus `0.000105` in coarse
regular cells. `off-first-step-localization.json` records the coarse/fine,
cut/regular breakdown. This control preceded the explicit `eb_noslip=0` input;
its implicit native no-slip switch is inactive because viscosity is disabled
(`Source/Diffusion.cpp` applies it only with `diffuse_vel`). The final
transport-on case explicitly sets matching slip conditions.

A debugger inspection of the first accepted step using the final input reports
`ok=true`, `theta=1`, and energy correction `2.225e-12 J`. It deliberately exits
before later steps and is not counted as another completed trajectory;
`first-step-limiter.log` retains the command result. This inspection gives no evidence that the overall limiter reduced the step.
The exact source of the discrepancy is not established. The next work should
compare the initial cut-face hydro flux/interpolation/divergence and coarse/fine
halo/reflux contributions before changing chemistry or widening any gate.

## Attempts and verification

All attempts remain under `build/tgreact-eb-amr-reference-20261019/runs/`:

- `matched_eb_amr_g2spehpn`: native startup refused the original 4³ grid;
  native EB requires blocking factor at least eight.
- `matched_eb_amr_s8s5yzui`: native stopped on NaNs after step one with
  zero-valued covered-body cells. The case now uses the native default
  computable body state, preserving the NaN rejection.
- `matched_eb_amr_tnfgo3i3`: native completed; the new Fortran adapter refused
  because its required held-energy flux argument was missing. The adapter
  was fixed without changing the production solver.
- `matched_eb_amr_i9x8uttw`: reused that native capture; Fortran completed but
  failed field gates. Native's implicit EB no-slip setting was then made
  explicitly slip to match the intended Fortran boundary.
- `matched_eb_amr_n7d04w93`: transport-off diagnostic completed, field gate failed.
- `matched_eb_amr_62p3fjsx`: final matching-slip transport-on pair completed,
  field gate failed as recorded above.

There are three completed Fortran trajectories, one refused adapter run,
three completed native trajectories (one reused), one native startup refusal
and one native NaN stop. There is no passing full comparison and no skipped
comparison recorded as passed. The debugger run is separate diagnostic work.

LiDryer provenance/regeneration, the new target build, Python syntax, project
contract and `git diff --check` pass. Unrelated unit suites/full matrices were
not repeated. Branch `main`, base
`0a0e00ebe9138145fbb0e0e745665d0d768df009`, with the pre-existing dirty tree
retained. GNU Fortran 13.3 Debug, Open MPI 4.1.6, CMake 3.28.3 and SUNDIALS 7.2;
no compiler/MPI/backend configuration changed. `source-identity.json`,
`source/`, `CMakeCache.txt`, build/static logs, run-specific commands and hashes
retain exact source and artifact identity. This does not qualify application
initialization, MPI, regridding, convergence or physical ignition.

The final source cleanup removes one trailing space from probe line 74.
`compared-source/` retains the exact preceding probe text used for the recorded
comparison; its hash and the different final rebuilt binary identity are
recorded in `verification.json`. The numerical comparisons were not repeated
for that whitespace-only rebuild.
