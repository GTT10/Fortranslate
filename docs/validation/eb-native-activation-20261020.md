# Native activation-basis compatibility, 2026-10-20

Matching the pinned native activation-energy conversion removes all six
longer-trajectory O2/H2O field failures in the held-EOS diagnostic. At 1.2
microseconds, their mass-fraction differences fall by about 118 and 115 times.
All field limits pass at all thirteen sampled times, with species/element
conservation and the independent numerical-energy budget intact.

This establishes the main cause of the accumulated chemical discrepancy in
this bounded case. It is **not** qualification against unchanged PeleC:
the native comparator still uses the separately identified held spatial EOS.
The canonical Cantera mechanism and its original failed comparisons remain
unchanged. Time/spatial convergence and physical validation remain open.

## Source-derived compatibility change

Pinned CEPTR `constants.py` defines the calorie gas constant as
`1.98721558317399615845 cal/mol/K`, which its `production.py` uses to produce
Arrhenius activation temperatures. This constant also appears as `ruc` in
the pinned LiDryer generated header. The canonical Cantera importer instead
converts cal/mol to J/kmol using 4184 and divides by its SI gas constant.
The resulting activation temperatures differ by a factor of
`0.9999943013061792`. This is a generator-convention difference, not a fitted
reaction multiplier.

`tools/prepare_lidryer_native_activation.py` checks the pinned YAML/header
SHA256 values against `references/pelec_tgreact_amr.json`, verifies the
canonical source and pinned Cantera environment, and checks all 17 nonzero
native activation coefficients against the source-derived conversion.
It generates a separate YAML with activation energies expressed in kelvin,
then uses the existing importer and generator. No generated Fortran rates
or original reference/source files are edited.

The bundle has 23 activation-energy entries, including both falloff limits
and zeros. Automatic post-generation comparison requires the source-derived
conversion for every entry and exact identity of all other reaction
parameters, species, thermo and transport data. Modified prefactors,
activation energies, stoichiometry and thermo were independently refused.
Earlier preparation checks also refused a changed native calorie constant,
activation coefficient and source SHA. Existing output directories are
refused by construction.

## Coupled result

The Fortran probe links the separately generated mechanism with the same
production EB/CVODE/runtime libraries and probe object used for the completed
time study. It runs twelve 100 ns root steps with recursive ratio-two MOL2,
independent initial reaction history and pinned Simple transport. The native
held-EOS snapshots from the completed longer run are reused byte for byte.
All 1024 initial gas/UEINT/T cells match exactly; comparison uses 544 composite
fluid cells. Geometry, transport, stepping controls and limits are unchanged.

| At 1.2 microseconds | Canonical Cantera basis | Native activation basis | Limit |
| --- | ---: | ---: | ---: |
| O2 mass-fraction difference | 5.419574964e-7 | 4.582725699e-9 | 1e-7 |
| H2O mass-fraction difference | 4.469964739e-7 | 3.896800586e-9 | 1e-7 |
| Temperature difference (K) | 2.299841849e-4 | 6.936792738e-5 | 1e-3 |
| Scaled total-energy difference | 3.420620610e-9 | 5.428088932e-11 | 1e-8 |
| Relative hydrogen-consumption difference | 9.846291550e-5 | 9.820404810e-7 | 1e-3 |

PeleF consumes 0.722680245% of the initial hydrogen, versus 0.722680955% in
the native diagnostic. Maximum PeleF mass/element drift is `5.9953e-14`,
and maximum relative energy-budget residual is `7.4163e-17`. Species remain
nonnegative with closed species sums. No tolerance was relaxed. The small
remaining discrepancies and the previously observed time order near 1.1--1.2
are not erased by this result.

## Reproduction and evidence

Artifacts are in `build/eb-native-activation-20261020/`. `run.json` records
the exact compile/link/run commands, object/library hashes, binary identity,
imported input and all reused native snapshot hashes. The completed run took
211.44 seconds with GNU Fortran 13.3, runtime/FPE checks and SUNDIALS 7.2.0.
No compiler, MPI or optional-backend configuration was changed; new mechanism
objects/modules were built in their own directory. `comparison.json` records
all thirteen times, zero diagnostic field exceedances and the energy budget.

`mechanism/activation-basis.json` records generation at run time.
`prepare_lidryer_native_activation_at_run.py` retains that exact tool version.
The subsequent guard extension was checked by regenerating in
`verified-regeneration/`: source YAML, bundle and Fortran bytes are identical
to those used by the completed run, so the CFD calculation was not repeated.
`generation-guard-checks.json` records the four additional expected refusals;
`import-checks.json` records the earlier three and the independent parameter
comparison. `verification.json` records final identities and completed checks.

The existing source-provenance/regeneration preflight, Python syntax,
`git diff --check` and project contract pass. One new coupled Fortran run
completed; the native trajectory was reused. No broad unit suite, MPI/restart
campaign, full matrix, installed-binary audit or physical-validation run was
repeated. The twelve failures against the original pinned spatial algorithm
remain separately recorded and open.
