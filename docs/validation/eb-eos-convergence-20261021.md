# Caloric inversion precision in the EB-AMR diagnostic, 2026-10-21

The tighter-chemistry 32-cubed prefix still amplified small z-momentum
differences through its shared coarse-to-fine limiter. Saved-state caloric
evaluation and an operator-only intervention isolated a contribution from
the pinned inverse's `1e-6 K` stopping rule; see the preceding
[forced-reactor and residual record](eb-forced-cvode-20261021.md).

## Completed private intervention

`build/eb-tight-eos-prefix-20261021/` completed five spatial evaluations
with private PeleF physical and signed-trial inverse stops changed to
`1e-12 K`. The native five-evaluation capture was reused unchanged with
its default inverse. Chemistry tolerances stayed at relative `1e-12`
and density absolute `1e-18 kg/m3`. Both executions used the corrected
forced-CVODE backend. This isolates the PeleF EOS precision intervention
from the preceding chemical-tolerance/backend changes.

| Maximum absolute difference | Default PeleF inverse | Tighter PeleF inverse |
| --- | ---: | ---: |
| First coarse predictor energy rate, J/m3/s | 669.872 | 11.2897 |
| First fine corrector halo energy, J/m3 | 8.92198 | 0.00975415 |
| First fine corrector halo x-momentum, kg/m2/s | 0.00831451 | 1.41381e-5 |
| First fine corrector energy rate, J/m3/s | 1.80497e6 | 3509.43 |
| PeleF capture wall time, s | 734.275 | 425.590 |

The z-momentum-controlled x limiter at coarse parent `(15,0,15)` changes
from 0.60627244 to 0.42469832, close to native 0.42500654. At `(13,14,18)`
it changes from 0.51055038 to 0.64059384, close to native 0.64058110.
Native caloric residuals remain: the second coarse corrector's temperature
gap is still `7.85e-7 K`, consistent with the preceding native residual audit.
The change therefore substantially reduces the diagnosed amplification but
does not prove exact parity or full-trajectory acceptance.

Wrapper exit 0, solver's intentional exit 86, completion sentinel, all five
hydro/transport captures and three comparison scripts were verified. The
same-thread completion notification is in `job/`; input/object hashes,
commands and elapsed times are in `capture-record.json`. The preceding
failed/timed-out captures are retained.

## Implementation and qualification boundary

`configure_nasa7_eos_policy` now accepts optional `converged_inverse=.true.`.
It selects the native `CONVERGENCE` settings: `1e-12 K` stopping tolerance
and 5000 iterations. Physical and bounded signed-trial inverses use this
mode consistently. The existing default remains `1e-6 K`, 200 iterations;
the library's separate strict inverse, fitted coefficients, boundary-Cv
extrapolation, signed composition bounds and endpoint rollback are unchanged.

The flag must be uniform across the mixture and requires pinned EOS mode.
Configuration is transactional. Its distinct `PELEF_FUEGO_EOS_CONVERGENCE_V1`
identity is used by the existing checkpoint and MPI policy comparisons;
default identity strings remain unchanged. No new MPI execution is claimed.

The EB reference probe accepts V4 with explicit density reaction tolerances
and a `0`/`1` EOS-mode row. Mesh diagnostic configuration can select
`fortran_eos_convergence: true`; reports distinguish PeleF and native modes,
and the mesh analyzer rejects mixed EOS precision. This is explicitly a
PeleF-only intervention against the already recorded native binary.

`build/eb-eos-convergence-20261021/` owns the clean Debug/CVODE build and
focused verification. The derived species type changed, so old objects are
not mixed into the new runtime. The clean build completes in 115.363 s with
GNU Fortran 13.3.0 (`/usr/bin/f95` alias), Debug bounds/FPE checks and SUNDIALS
7.2. Its four focused units pass in 1.29 s: analytic physical and signed
inverses, forced-CVODE admission/rollback and regular checkpoint precision
mismatch/same-mode continuation. These are four executables, not a full suite.

The default 8-cubed probe produces six CSV files byte-identical to the saved
pre-change fixed-backend run. The opt-in 8-cubed coupled comparison passes all
unchanged field, inventory and budget limits at 0, 50 and 100 ns with explicit
relative tolerance `1e-12` and density absolute tolerance `1e-18 kg/m3`.
Native takes 17.475 s and PeleF 41.404 s at these tighter chemical tolerances.
The completed report is
`mesh8/matched_eb_amr_apqmwnxt/checks.json`, with both EOS modes recorded.

The two existing mesh-analysis Python checks pass, including refusal of mixed
EOS precision. Three malformed study controls and three invalid V4 mode rows
are refused. Canonical/variant provenance and regeneration, Python syntax,
`git diff --check` and the project contract pass. Exact commands/logs are in
`configure.json`, `build.json`, `verification.json`, `control-refusals.json`
and `source-checks.json`; `source-identity.json` records current main base
`0a0e00ebe9138145fbb0e0e745665d0d768df009`, toolchain and source hashes.
The first verification driver assumed a `gfortran` command name and stopped
at the link parser; it was corrected to read CMake's `f95` alias, preserving
and reusing the already completed checks. That partial record is retained as
`verification-at-link-parser-failure.json`.

## Completed follow-up

`build/eb-eos-mesh-study-20261021/` reuses the completed 8-cubed report and
runs matching 16- and 32-cubed trajectories through 100 ns at common 50-ns
steps. Both complete within their 1800-second solver limits and pass all
unchanged comparison gates. The three passing reports enter the conservative
mesh analyzer; the campaign finishes in 42.82 min. Its observed orders and the
subsequent initial-profile diagnosis are in the
[three-grid record](eb-three-grid-20261021.md).

This configured 32-cubed field comparison is complete. Asymptotic convergence,
unchanged native parity, the configuration/install matrix, MPI/restart
campaigns and physical validation remain open.
